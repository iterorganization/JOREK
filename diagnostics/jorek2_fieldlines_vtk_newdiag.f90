!**************************************************************************************************************************
!* This diagnostics evaluates arbitrary physical expressions (using the new_diag framework)                               *
!* along field lines and writes this information to both vtk and txt files.                                               *
!* It uses the file jorek_restart.h5 to load the simulation state,                                                        *
!* and a text file called stpts containing required information                                                           *
!* about the field lines along which physical expressions are to be to evaluated,                                         *
!* and the physical expressions to be evaluated.                                                                          *
!* Output is generated in the binary file field_lines.vtk (using a coordinate system consistent with jorek2vtk)           *
!* and ASCII file field_lines.txt (using the coordinate system in https://www.jorek.eu/wiki/doku.php?id=coordinates))     *
!*                                                                                                                        *
!* For details, see:                                                                                                      *
!* https://www.jorek.eu/wiki/doku.php?id=jorek2_fieldlines_vtk_newdiag                                                    *
!**************************************************************************************************************************

program jorek2_fieldlines_vtk_newdiag

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
use mod_new_diag
use mod_impurity, only : init_imp_adas

implicit none

real*8,allocatable  :: rp(:), zp(:), R_all(:), Z_all(:), C_all(:), R_strike(:), Z_strike(:), P_strike(:), C_strike(:)
real*4,allocatable  :: Xfield(:,:), Yfield(:,:), Zfield(:,:), Pfield(:,:), Tfield(:,:,:)
real*4,allocatable  :: XfieldVTK(:,:), YfieldVTK(:,:), ZfieldVTK(:,:)
integer,allocatable :: Nfield(:)
character*12, allocatable :: scalar_names(:), vector_names(:)
integer :: i, j, iside_i, iside_j, ip, i_line, n_lines, i_tor, i_harm, i_var_psi, i_dir, k, m, ntotal
integer :: n_start, n_end, n_large, n_total, n_delta, n_scalars, n_vectors
integer :: i_elm, ifail, i_phi, n_phi, i_turn, i_elm_out, i_elm_prev, i_elm_tmp,i_steps, i_field
real*8  :: R_line, Z_line, s_line, t_line, p_line, s_mid, t_mid, p_mid, s_out, t_out
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: tol, delta_phi, Zjac, psi_s, psi_t, R_in, Z_in, R_out, Z_out, Rmin, Rmax, Zmin, Zmax, delta_s, delta_t, R_keep, Z_keep
real*8  :: small_delta, small_delta_s, small_delta_t, delta_phi_local, delta_phi_step, total_phi
real*8  :: Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t, dl2, total_length, length_max, s_ini, t_ini, value_out
real*8  :: psi_norm_out, theta_out

character :: buffer*80, lf*1, str1*12, str2*12, str3*24
integer :: ivtk, i_var, my_id, ierr
logical :: psi_theta
real*8  :: coord_min(2), coord_max(2), coord_out(2)

character(len=512) :: s
character(len=12) :: name, expr
integer :: nr, ntour, curr, n_dir
real*8  :: rr, zz, psi
type(t_pol_pos_list) :: pol_pos_list
type(t_tor_pos_list) :: tor_pos_list
type(t_expr_list)    :: expr_list
integer, allocatable :: n_turn(:)
real*8, allocatable :: result(:,:,:,:), res0d(:)
real*8, allocatable :: R_start(:), Z_start(:), P_start(:)
character*12, allocatable :: scalar_exprs(:)

write(*,*) '***************************************'
write(*,*) '* JOREK2_fieldlines_vtk_newdiag       *'
write(*,*) '***************************************'
write(*,*) ' nperiod : ',n_period

my_id=0

call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Read start points and expressions from file 'stpts'.
!
! Example for a stpts file:
!   +-------------------------------------------------------
!   |# n_scalars
!   |   3
!   |# expr  name
!   |  Te    Te_eV
!   |  ne    ne_m-3
!   |  nimp  nimp_m-3
!   |# n_lines
!   |   100
!   |# nr    R_start   Z_start   phi_start   n_turns
!   |   1    3.003     1.606     6.0900      1
!   |   2    2.988     1.558     6.0900      1
!   |   3    2.990     1.560     6.0900      1
!   |   4    2.946     1.576     6.0900      1
!   |   5    2.989     1.571     6.0900      1
!   |   6    2.995     1.649     6.0900      1
!   |   7    3.002     1.648     6.0900      1
!   |   8    2.977     1.550     6.0900      1
!   |   9    3.017     1.631     6.0900      1
!   |  10    3.033     1.621     6.0900      1
!   | ...
!   | 100    2.997     1.620     6.0900      1
!   +-------------------------------------------------------
!
! n_scalars is the number of scalar expressions to be evaluated,
! expr and name are the expressions (as defined in the new_diag framework) and its name,
! n_lines are the number of field lines to be evaluated
! (in both backward and forward direction so the final number of computed field lines will be twice this number).
! nr is the field line number, R_start, Z_start and phi_start its initial coordinate,
! and n_turns the number of toroidal turns along which the field line should be followed.
! As in the jorek2_poincare diagnostics, the user can either provide the initial coordinates of all the field lines,
! or specify groups of field lines where only the first and last lines are explicitly given,
! and their initial coordinates are assumed to be uniformily distributed between the first and the last.
!
open(21, file='stpts', status='old', action='read', iostat=ierr)

if ( ierr == 0 ) then ! stpts file exists, use it.

  read(21, '(a)') s ! read comment line (ignored)
  read(21,*) n_scalars
  if ( n_scalars < 1 ) then
    write(*,*) 'ERROR in stpts file: n_scalars must be >= 1.'
    stop
  end if

  allocate(scalar_names(n_scalars),scalar_exprs(n_scalars))

  read(21, '(a)') s ! read comment line (ignored)
  do i = 1, n_scalars
     read(21,*) expr,name
     scalar_exprs(i)=expr
     scalar_names(i)=name
  end do

  write (*,*) 'Scalar expressions: ',scalar_exprs
  write (*,*) 'Scalar names:       ',scalar_names
  
  read(21, '(a)') s ! read comment line (ignored)
  read(21,*) n_lines
  if ( n_lines < 1 ) then
    write(*,*) 'ERROR in stpts file: n_lines must be >= 1.'
    stop
  end if

  write (*,*) n_lines
  
  read(21, '(a)') s ! read comment line (ignored)

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
      
      do i_line = curr + 1, nr
        R_start(i_line) = R_start(curr) + ( rr - R_start(curr) ) * ( real(i_line-curr) / real(nr-curr) )
        Z_start(i_line) = Z_start(curr) + ( zz - Z_start(curr) ) * ( real(i_line-curr) / real(nr-curr) )
        P_start(i_line) = P_start(curr) + ( psi- P_start(curr) ) * ( real(i_line-curr) / real(nr-curr) )
        n_turn(i_line)  = nint( n_turn(curr) + real( ntour - n_turn(curr) ) * ( real(i_line-curr) / real(nr-curr) ) )
      end do
      
    end if
    
    curr = nr
    
  end do
  
  close(21)

else ! if no stpts file exists, use the following hard-coded default startpoints

  n_scalars=3
  scalar_exprs = (/'Te  ','ne  ','nimp'/)
  scalar_names = (/'Te_eV   ','ne_m-3  ','nimp_m-3'/)

  n_lines = 50

  allocate( R_start(n_lines), Z_start(n_lines), P_start(n_lines), n_turn(n_lines) )

  do i_line = 1, n_lines
    R_start(i_line) = 1.7156 + (2.18-1.7156) * float(i_line-1)/float(n_lines-1)
    Z_start(i_line) = 0.12237
    P_start(i_line) = 0.d0
  end do

  n_turn  = 500
 
end if

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  if (my_id .eq. 0 ) then
     write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
  endif
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                                       ! define the basis functions at the Gaussian points

! --- Initialize the new_diag framework and print some information (.true.)
! additional to jorek2_fieldlines_vtk
call init_new_diag(.true.)
expr_list = exprs(scalar_exprs, n_scalars)

! --- Read ADAS data and generate coronal equilibrium is needed
! additional to jorek2_fieldlines_vtk
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  call init_imp_adas(0)
#endif


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
n_phi   = 1000
n_dir   = 2 ! forward and backward directions
!n_dir   = 1 ! only forward direction

delta_phi = 2.d0 * PI / float(n_period*n_phi)
tol       = 1.d-6

i_var_psi = 1

n_start = 1
n_end   = element_list%n_elements

n_large = maxval(n_turn) * n_phi * 10

write(*,*) ' n_lines : ',n_lines,n_large

allocate(R_strike(n_dir*n_lines),Z_strike(n_dir*n_lines),P_strike(n_dir*n_lines),C_strike(n_dir*n_lines))

allocate(R_all(n_dir*n_lines),Z_all(n_dir*n_lines),C_all(n_dir*n_lines))

allocate(Xfield(n_large,n_dir*n_lines),Yfield(n_large,n_dir*n_lines),Zfield(n_large,n_dir*n_lines),Pfield(n_large,n_dir*n_lines))
allocate(XfieldVTK(n_large,n_dir*n_lines),YfieldVTK(n_large,n_dir*n_lines),ZfieldVTK(n_large,n_dir*n_lines))
allocate(Nfield(n_dir*n_lines),Tfield(n_large,n_dir*n_lines,n_scalars))

R_all    = 0.d0; Z_all    = 0.d0; C_all = 0.d0
R_strike = 0.d0; Z_strike = 0.d0; P_strike = 0.d0; C_strike = 0.d0
Xfield   = 0.d0; Yfield   = 0.d0; Zfield = 0.d0; Pfield = 0.d0; Tfield = 0.d0; Nfield = 0
XfieldVTK   = 0.d0; YfieldVTK   = 0.d0; ZfieldVTK = 0.d0

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

do i=1,n_lines

   write(*,*)
   write(*,'(1x,2(a,i6),a,2f8.3)') 'Line',i,' of',n_lines,' started at',R_start(i),Z_start(i)

  do i_dir = -1*(-1)**n_dir,1,n_dir

     write(*,'(1x,a,i6)') 'Direction',i_dir

     call find_RZ(node_list,element_list,R_start(i),Z_start(i),R_out,Z_out,i_elm,s_out,t_out,ifail)
 
     if (ifail .ne. 0) write(*,*) "Can not find RZ,", ifail 
     if (ifail .ne. 0) exit

     R_line = R_start(i)
     Z_line = Z_start(i)
     p_line = P_start(i)
     s_line = s_out
     t_line = t_out

      i_field = 0
      i_line  = i_line + 1

      delta_phi = 2.d0 * PI * float(i_dir) / float(n_period*n_phi)

      total_length = 0.d0
      total_phi    = p_line

      Xfield(1,i_line) = R_line * cos(total_phi)
      Yfield(1,i_line) = -R_line * sin(total_phi)
      Zfield(1,i_line) = Z_line
      XfieldVTK(1,i_line) = R_line * cos(total_phi-P_start(i))
      ZfieldVTK(1,i_line) = R_line * sin(total_phi-P_start(i))
      YfieldVTK(1,i_line) = Z_line
      Pfield(1,i_line) = total_phi

      call create_pol_pos(pol_pos_list, ierr, node_list, element_list, ES, ielm=i_elm, s=s_line, t=t_line)
      call create_tor_pos(tor_pos_list, ierr, phi=total_phi)
      call eval_expr(ES, SI_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
      call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)
      Tfield(1,i_line,:) = res0d
      
      i_field = 1

      R_all(i_line) = R_line
      Z_all(i_line) = Z_line

      do i_turn = 1, n_turn(i)
         if ( mod(i_turn-1,max(n_turn(i)/6+1,5)) == 0 ) then
            write(*,'(3x,2(a,i6))') 'Turn',i_turn,' of',n_turn(i)
         end if
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
	    Yfield(i_field,i_line) = -R_in * sin(total_phi)
	    Zfield(i_field,i_line) = Z_in
	    XfieldVTK(i_field,i_line) = R_in * cos(total_phi-P_start(i))
	    ZfieldVTK(i_field,i_line) = R_in * sin(total_phi-P_start(i))
	    YfieldVTK(i_field,i_line) = Z_in
	    Pfield(i_field,i_line) = total_phi
	    Nfield(i_line)         = Nfield(i_line) + 1

            call create_pol_pos(pol_pos_list, ierr, node_list, element_list, ES, ielm=i_elm, s=s_line, t=t_line)
            call create_tor_pos(tor_pos_list, ierr, phi=total_phi)
            call eval_expr(ES, SI_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr)
            call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)
            Tfield(i_field,i_line,:) = res0d

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

enddo ! end of main loop

n_lines = i_line

!open(20,file='fieldlines.txt')
!do i_line=1,n_lines
!   write(20,'(i8,8e16.8)') i_line,C_all(i_line),R_all(i_line),Z_all(i_line)
!enddo
!close(20)

open(20,file='field_lines.txt')
write(20,'(13A16)') 'line', 'x', 'y', 'z', 'Phi', (scalar_names(i_var),i_var=1,n_scalars)
do i=1,n_lines
   do j=1,Nfield(i)
      write(20,'(i16,12e16.8)') i,Xfield(j,i),Yfield(j,i),Zfield(j,i),Pfield(j,i),(Tfield(j,i,i_var),i_var=1,n_scalars)
   enddo
enddo
close(20)

!do i_line=1,n_lines
!   write(*,'(i8,8e12.4)') i_line,C_strike(i_line),R_strike(i_line),Z_strike(i_line),P_strike(i_line)
!enddo

!--------------------------------------------------- write the binary VTK file
!n_scalars=1
!### call tr_allocate(scalar_names,1,n_scalars,"scalar_names")
!allocate(scalar_names(1))
!scalar_names = (/ 'T           '/)

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

write(ivtk) (([XfieldVTK(j,i),YfieldVTK(j,i),ZfieldVTK(j,i)], j=1,Nfield(i)),i=1,n_lines)

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
  write(ivtk) (((Tfield(j,i,i_var)), j=1,Nfield(i)),i=1,n_lines)
enddo

!do i_var =1, n_vectors
!  buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf                      ; write(ivtk) trim(buffer)
!  write(ivtk) ((vectors(j,i,i_var),i=1,3),j=1,nnos)
!enddo

close(ivtk)

end program jorek2_fieldlines_vtk_newdiag



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
