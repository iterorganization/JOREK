program jorek2_fieldlines_vtk

use constants
use data_structure
use phys_module
use basis_at_gaussian                                       
use elements_nodes_neighbours
use tr_module
use mod_neighbours
use mod_import_restart
use equil_info, only : get_psi_n, ES
use mod_interp

implicit none

real*8,allocatable  :: rp(:), zp(:), R_all(:), Z_all(:), C_all(:), R_strike(:), Z_strike(:), P_strike(:), C_strike(:)
real*4,allocatable  :: Xfield(:,:), Yfield(:,:), Zfield(:,:), Tfield(:,:)
integer,allocatable :: Nfield(:)
character*12, allocatable :: scalar_names(:), vector_names(:)
integer :: i, j, iside_i, iside_j, ip, i_line, n_lines, i_tor, i_harm, i_var_psi, i_dir, k, m, nr, np, ntotal
integer :: n_start, n_end, n_large, n_total, n_delta, n_scalars, n_vectors
integer :: i_elm, ifail, i_phi, n_phi, i_turn, n_turns, i_elm_out, i_elm_prev, i_elm_tmp,i_steps, i_field
real*8  :: R_start, Z_start, P_start, R_line, Z_line, s_line, t_line, p_line, s_mid, t_mid, p_mid, s_out, t_out
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: tol, delta_phi, Zjac, psi_s, psi_t, R_in, Z_in, R_out, Z_out, Rmin, Rmax, Zmin, Zmax, delta_s, delta_t, R_keep, Z_keep
real*8  :: small_delta, small_delta_s, small_delta_t, delta_phi_local, delta_phi_step, total_phi
real*8  :: Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t, dl2, total_length, length_max, s_ini, t_ini, value_out
real*8  :: psi_norm_out, theta_out

character :: buffer*80, lf*1, str1*12, str2*12, str3*24
integer :: ivtk, i_var, my_id, ierr
logical :: psi_theta
real*8  :: coord_min(2), coord_max(2), coord_out(2)

namelist /fieldlines_vtk_params/ psi_theta, n_turns, n_phi, n_lines, coord_min, coord_max

write(*,*) '***************************************'
write(*,*) '* JOREK2_fieldlines_vtk               *'
write(*,*) '***************************************'
write(*,*) ' nperiod : ',n_period

my_id=0

call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters 
  ! steps
n_turns = 20             ! number of toroidal turns to follow a fieldline
n_phi   = 1000           ! number of steps per toroidal turn
n_lines = 800
  !default coordinates of the starting points.
psi_theta = .true.
coord_min(1)= 0.95
coord_min(2)= -PI/20.d0
coord_max(1)= 0.97
coord_max(2)= +PI/20.d0

! --- Read parameters from namelist file 'fieldlines_vtk.nml' if it exists
open(42, file='fieldlines_vtk.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
if (my_id .eq. 0 ) then
   write(*,*) 'Reading parameters from fieldlines_vtk.nml namelist.'
endif
read(42,fieldlines_vtk_params)
close(42)
end if

if (my_id .eq. 0 ) then
   write(*,*)
   write(*,*) 'Parameters:'
   write(*,*) '-----------'
   write(*,*) 'n_turns = ', n_turns
   write(*,*) 'n_phi = ', n_phi
   write(*,*) 'n_lines = ', n_lines
   write(*,*) 'psi_theta = ', psi_theta
   write(*,*) 'coord_min = ', coord_min
   write(*,*) 'coord_max = ', coord_max
endif


do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  if (my_id .eq. 0 ) then
     write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
  endif
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                                       ! define the basis functions at the Gaussian points

allocate(element_neighbours(4,element_list%n_elements))

element_neighbours = 0

do i=1,element_list%n_elements

  do j=i+1,element_list%n_elements

    if (neighbours(node_list, element_list%element(i),element_list%element(j),iside_i,iside_j)) then
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

!n_turns = 20
!n_phi   = 1000

np = 5
nr = 2

delta_phi = 2.d0 * PI / float(n_period*n_phi)
tol       = 1.d-6

i_var_psi = 1

n_start = 1
n_end   = element_list%n_elements

!n_lines = 800 ! 2* (n_end - n_start + 1) * nr * np
n_large = n_turns * n_phi * 10

write(*,*) ' n_lines : ',n_lines,n_large

allocate(R_strike(n_lines),Z_strike(n_lines),P_strike(n_lines),C_strike(n_lines))

allocate(R_all(n_lines),Z_all(n_lines),C_all(n_lines))

allocate(Xfield(n_large,n_lines),Yfield(n_large,n_lines),Zfield(n_large,n_lines),Nfield(n_lines),Tfield(n_large,n_lines))

R_all    = 0.d0; Z_all    = 0.d0; C_all = 0.d0
R_strike = 0.d0; Z_strike = 0.d0; P_strike = 0.d0; C_strike = 0.d0
Xfield   = 0.d0; Yfield   = 0.d0; Zfield = 0.d0;  Tfield = 0.d0; Nfield = 0

Rmin = 1.d20; Rmax = -1.d20; Zmin = 1.d20; Zmax=-1.d20
do i=1,node_list%n_nodes
  Rmin = min(Rmin,node_list%node(i)%x(1,1,1))
  Rmax = max(Rmax,node_list%node(i)%x(1,1,1))
  Zmin = min(Zmin,node_list%node(i)%x(1,1,2))
  Zmax = max(Zmax,node_list%node(i)%x(1,1,2))
enddo

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

i_line = 0

write(*,*) ' number of elements : ',element_list%n_elements

if (psi_theta) then
   !------------------------------------------------- find x-point(s), psi_bnd and psi_axis
   if (my_id .eq. 0 ) then
      write(*,*) ' xcase,1st x-point:R,Z,psi: ',xcase, ES%R_xpoint(1),ES%Z_xpoint(1),ES%psi_xpoint(1),ES%psi_bnd
      !   write(*,*) ' PSI_XPOINT : ',psi_xpoint,i_elm_xpoint
      write(*,*) ' PSI_AXIS : ', ES%psi_axis,ES%i_elm_axis
      write(*,*) ' RZ_AXIS : ',  ES%R_axis, ES%Z_axis
   endif
endif

do i =n_start, n_end

  if ( i_line .ge. n_lines ) cycle

  do k=1, nr

    if ( i_line .ge. n_lines ) cycle
    
    s_ini = real(k)/real(nr+1)

    do m=1, np

      if ( i_line .ge. n_lines ) cycle
            
      t_ini = real(m)/real(np+1)

      call var_value(i,1,s_ini,t_ini,0.d0,value_out)

      call interp_RZ(node_list,element_list,i,s_ini,t_ini,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
      psi_norm_out = get_psi_n(value_out, Z_out)
      theta_out = atan2( (Z_out - ES%Z_axis) , (R_out - ES%R_axis) ) !/ (2.d0*PI)      
      if (psi_theta) then
         coord_out(1) = psi_norm_out
         coord_out(2) = theta_out
      else
         coord_out(1) = R_out
         coord_out(2) = Z_out
      endif
         if ((coord_out(1) .gt. coord_min(1)) .and. (coord_out(1) .lt. coord_max(1)) .and. &
              (coord_out(2) .gt. coord_min(2)) .and. (coord_out(2) .lt. coord_max(2)) .and. (i_line .lt. n_lines) ) then
	  
	  write(*,'(i6,2f7.3)') i_line, coord_out

      do i_dir = -1,1,2

      i_field = 0
      i_line  = i_line + 1

      s_line = s_ini
      t_line = t_ini

      delta_phi = 2.d0 * PI * float(i_dir) / float(n_period*n_phi)

      !call interp_RZ(node_list,element_list,i,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
             

      total_length = 0.d0
      total_phi    = 0.d0

      i_elm = i
      R_start = R_out
      Z_start = Z_out
      P_start = 0.d0

      Xfield(1,i_line) = R_start * cos(total_phi)
      Zfield(1,i_line) = R_start * sin(total_phi)
      Yfield(1,i_line) = Z_start

      call var_value(i_elm,6,s_line,t_line,total_phi,value_out)
      Tfield(1,i_line) = value_out
      
      i_field = 1

      R_all(i_line) = R_start
      Z_all(i_line) = Z_start

      R_line = R_start
      Z_line = Z_start
      p_line = P_start

      do i_turn = 1, n_turns

        do i_phi=1,n_phi

          delta_phi_local = 0.d0

          i_steps = 0

          do while ((abs(delta_phi_local) .lt. abs(delta_phi)) .and. (i_steps .lt.10) )

            i_steps = i_steps + 1

            delta_phi_step = delta_phi - delta_phi_local

!	write(*,'(5i6,3e16.8)') i_line,i_turn,i_phi,i_steps,i_elm,delta_phi_step, delta_phi, delta_phi_local

            call step(i_elm,s_line,t_line,p_line,delta_phi_step,delta_s,delta_t,R,Z,R_s,R_t,Z_s,Z_t)

            s_mid = s_line + 0.5d0 * delta_s
            t_mid = t_line + 0.5d0 * delta_t
            p_mid = p_line + 0.5d0 * delta_phi_step

            call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t,Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t)

            small_delta_s = 1.d0

            if  (s_line + delta_s .gt. 1.d0) then

              small_delta_s = (1.d0 - s_line)/delta_s

!	      write(*,*) ' 1 : ',small_delta_s,s_line,delta_s

            elseif  (s_line + delta_s .lt. 0.d0) then

              small_delta_s = abs(s_line/delta_s)

!	      write(*,*) ' 2 : ',small_delta_s,s_line,delta_s

            endif

            small_delta_t = 1.d0

            if  (t_line + delta_t .gt. 1.d0)  then

              small_delta_t = (1.d0 - t_line)/delta_t

!	      write(*,*) ' 3 : ',small_delta_t,t_line,delta_t

            elseif  (t_line + delta_t .lt. 0.d0)  then

              small_delta_t = abs(t_line/delta_t)

!	      write(*,*) ' 4 : ',small_delta_t,t_line,delta_t

            endif

            small_delta = min(small_delta_s, small_delta_t)

!           write(*,'(A,5e16.8)') ' small delta : ',small_delta,delta_s,delta_t,s_line,t_line

            if (small_delta .lt. 1.d0)  then

              s_mid = s_line + 0.5d0 * small_delta * delta_s
              t_mid = t_line + 0.5d0 * small_delta * delta_t
              p_mid = p_line + 0.5d0 * small_delta * delta_phi_step

             call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t,Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t)

             if (small_delta_s .lt. small_delta_t) then

               if (s_line + delta_s .gt. 1.d0) then

	          s_line = 1.d0
                  t_line = t_line + small_delta * delta_t
                  p_line = p_line + small_delta * delta_phi_step

	          dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	          dl2 = dl2 * small_delta

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	          i_elm_prev = i_elm
                  i_elm      = element_neighbours(2,i_elm_prev)

	          if (i_elm .ne. 0) then

                    i_elm_tmp  = element_neighbours(4,i_elm)

	            if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (1)'

                    s_line = 0.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	            if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	              write(*,'(A,2i6,4f8.4)') ' error in element change (1) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

	          endif

                elseif (s_line + delta_s .lt. 0.d0) then

                  s_line = 0.d0
	          t_line = t_line + small_delta * delta_t
  	          p_line = p_line + small_delta * delta_phi_step

	          dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	          dl2 = dl2 * small_delta

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	          i_elm_prev = i_elm

	          if (i_elm .ne. 0) then

                    i_elm      = element_neighbours(4,i_elm_prev)
                    i_elm_tmp  = element_neighbours(2,i_elm)

                    if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (2)'

                    s_line = 1.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	            if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
	              write(*,'(A,2i6,4f8.4)') ' error in element change (2) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

	          endif

	        endif

	      else

                if (t_line + delta_t .gt. 1.d0) then

                  s_line = s_line + small_delta * delta_s
                  t_line = 1.d0
                  p_line = p_line + small_delta * delta_phi_step

	          dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	          dl2 = dl2 * small_delta

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                  i_elm_prev = i_elm
                  i_elm      = element_neighbours(3,i_elm_prev)

	          if (i_elm .ne. 0) then

                    i_elm_tmp  = element_neighbours(1,i_elm)

	            if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (3)'

                    t_line = 0.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	            if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
	              write(*,'(A,2i6,4f8.4)') ' error in element change (3) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

	          endif

                elseif (t_line + delta_t .lt. 0.d0) then

                  s_line = s_line + small_delta * delta_s
	          t_line = 0.d0
	          p_line = p_line + small_delta * delta_phi_step

	          dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
	          dl2 = dl2 * small_delta

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	          i_elm_prev = i_elm
                  i_elm      = element_neighbours(1,i_elm_prev)

	          if (i_elm .ne. 0) then

                    i_elm_tmp  = element_neighbours(3,i_elm)

	            if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (4)'

                    t_line = 1.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

	            if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
	              write(*,'(A,2i6,4f8.4)') ' error in element change (4) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

	          endif

                endif

	      endif

            else

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

	    i_field = i_field + 1
	    Xfield(i_field,i_line) = R_in * cos(total_phi)
	    Zfield(i_field,i_line) = R_in * sin(total_phi)
	    Yfield(i_field,i_line) = Z_in
	    Nfield(i_line)         = Nfield(i_line) + 1

            call var_value(i_elm,6,s_line,t_line,total_phi,value_out)
	    
	    Tfield(i_field,i_line) = value_out

            if (i_elm .eq. 0) exit

          enddo

          if (i_elm .eq. 0) exit

          if (i_steps .gt. 8) write(*,'(A,5i6)') ' WARNING : isteps ',i_line,i_turn,i_phi,i_steps,i_elm

        enddo ! end of a 2Pi turn (or before if end of open field line)

        if (i_elm .eq. 0) exit


      enddo  ! end of loop over toroidal turns

      if (i_dir .eq. -1) then
        C_all(i_line) = total_length
      else
        C_all(i_line) = min(C_all(i_line),total_length)
      endif

      enddo  ! end of two directions

!      if (i_elm .eq. 0) then
!        R_strike(i_line) = R_line
!        Z_strike(i_line) = Z_line
!        P_strike(i_line) = p_line
!        C_strike(i_line) = total_length
!      endif

    endif ! area selection

    enddo
  enddo


enddo ! end of loop over elements


n_lines = i_line

open(20,file='fieldlines.txt')
do i_line=1,n_lines
   write(20,'(i8,8e16.8)') i_line,C_all(i_line),R_all(i_line),Z_all(i_line)
enddo
close(20)

!do i_line=1,n_lines
!   write(*,'(i8,8e12.4)') i_line,C_strike(i_line),R_strike(i_line),Z_strike(i_line),P_strike(i_line)
!enddo

!--------------------------------------------------- write the binary VTK file
n_scalars=1
!### call tr_allocate(scalar_names,1,n_scalars,"scalar_names")
allocate(scalar_names(1))
scalar_names = (/ 'T           '/)


n_total = 0
do i=1,n_lines
  n_total = n_total + Nfield(i)
enddo
write(*,*) ' n_total = ',n_total

lf = char(10) ! line feed character

ivtk = 20

#ifdef IBM_MACHINE
open(unit=ivtk,file='field_lines.vtk',form='unformatted',access='stream',status='replace')
#else
open(unit=ivtk,file='field_lines.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN',status='replace')
#endif

buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
buffer = 'DATASET POLYDATA'//lf//lf                                                   ; write(ivtk) trim(buffer)

! POINTS SECTION
write(str1(1:12),'(i12)') n_total
buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)

write(ivtk) (([Xfield(j,i),Yfield(j,i),Zfield(j,i)], j=1,Nfield(i)),i=1,n_lines)

!LINES SECTION
write(str3(1:24),'(2i12)') n_lines,n_total+n_lines
buffer = 'LINES '//str3//lf                                                           ; write(ivtk) trim(buffer)

n_start = 0
do i=1,n_lines
  write(ivtk) Nfield(i),(n_start+j-1,j=1,Nfield(i))
  n_start = n_start + Nfield(i)
enddo

! POINT_DATA SECTION
write(str1(1:12),'(i12)') n_total
buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)

do i_var =1, n_scalars
  buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
  write(ivtk) (((Tfield(j,i)), j=1,Nfield(i)),i=1,n_lines)
enddo

!do i_var =1, n_vectors
!  buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf                      ; write(ivtk) trim(buffer)
!  write(ivtk) ((vectors(j,i,i_var),i=1,3),j=1,nnos)
!enddo

close(ivtk)

end program jorek2_fieldlines_vtk



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
use mod_interp

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
