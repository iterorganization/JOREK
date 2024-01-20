!**********************************************************************
!* program to convert a JOREK2 restart file into binary VTK format    *
!**********************************************************************

program jorek_diagnostics
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use mod_parameters, only: n_var
use data_structure
use phys_module
use basis_at_gaussian
use diffusivities, only: get_dperp, get_zkperp
use nodes_elements
use mod_boundary
use mod_import_restart
use equil_info, only : get_psi_n, ES
use mod_interp

implicit none

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors, my_id
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
integer,allocatable   :: ien (:,:)
integer               :: i, j, k, m, etype, ivtk, irst, int, i_var, i_tor, index, index_node, n_points
integer               :: n_bnd, iv, iv1, iv2, inode1, inode2, ierr, i_elm_xpoint(2), ifail, i_elm_axis, inode_avg
character             :: buffer*80, lf*1, str1*8, str2*8
character*8, allocatable :: scalar_names(:), vector_names(:)
real*4                :: float
real*8                :: sg, tg, phi, angle
real*8,allocatable    :: avg6(:), avg7(:), avg8(:), avg11(:), Ravg(:), Zavg(:)
real*8,allocatable    :: prf1(:), prf2(:), prf3(:), prf4(:), prf5(:), prf6(:), prf7(:), prf8(:), prf11(:), Rprf(:), Zprf(:)
logical,allocatable   :: tobedone(:)
real*8                :: BigR, xjac, psi_norm
real*8                :: R,R_s,R_t,R_st,R_ss,R_tt, Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8                :: PS,PS_s,PS_t,PS_st,PS_ss,PS_tt, VP,VP_s,VP_t,VP_st,VP_ss,VP_tt
real*8                :: RH,RH_s,RH_t,RH_st,RH_ss,RH_tt, TT,TT_s,TT_t,TT_st,TT_ss,TT_tt
real*8                :: rho, rho_s, rho_t, T, T_s, T_t, T_p, vpar, D_prof, ZK_prof, ZKpar_T
real*8                :: psi, psi_s, psi_t, psi_x, psi_y, T_x, T_y, BB2, normal, rho_norm, t_norm
real*8                :: distance, R_start, Z_start
logical               :: periodic
logical               :: without_n0_mode
integer               :: i_bnd_node, i_node

namelist /vtk_params/ nsub, periodic, without_n0_mode

write(*,*) '*********************************'
write(*,*) '* jorek_target2vtk              *'
write(*,*) '*********************************'

! --- Initialise input parameters and read the input namelist.
my_id = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub      = 5             ! Number of subdivisions of the cubic finite elements into linear pieces
periodic  = .false.
without_n0_mode = .false. ! If .true., the axisymmetric part will not be included

! --- Read parameters from namelist file 'vtk.nml' if it exists
open(42, file='vtk.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from vtk.nml namelist.'
  read(42,vtk_params)
  close(42)
end if
write(*,*)
write(*,*) 'Parameters:'
write(*,*) '-----------'
write(*,*) 'nsub            =', nsub
write(*,*) 'periodic        =', periodic
write(*,*) 'without_n0_mode =', without_n0_mode
write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*) 'F0        =', F0
write(*,*)
call flush_it(6)

ivtk = 21                 ! an arbitrary unit number for the VTK output file

n_scalars = 13             ! number of scalars to write to the VTK output file
n_vectors = 3

allocate(scalar_names(n_scalars), vector_names(n_vectors))

scalar_names = (/ 'flux    ','density ','T       ','Vpar    ','nV.n    ','nTV.n   ','KparT.n ', &
                  'Kperp.T ','Dperp.n ','B.n     ','nvT_gam ','n       ','nV3.n   '/)

vector_names = (/ 'B_field ','Velocity','normal  '/)

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                              ! define the basis functions at the Gaussian points

rho_norm = central_density*1.d20 * central_mass * 1.67d-27
t_norm   = sqrt(MU_zero*rho_norm)

! --- Find the lowest point on the outer divertor target (required later)
call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
R_start = 1.d99
Z_start = 1.d99
do i_bnd_node = 1, bnd_node_list%n_bnd_nodes
  i_node = bnd_node_list%bnd_node(i_bnd_node)%index_jorek
  R = node_list%node(i_node)%x(1,1,1)
  Z = node_list%node(i_node)%x(1,1,2)
  if ( ( node_list%node(i_node)%boundary == 3 ) .and. ( R > ES%R_xpoint(1) ) .and. ( Z < Z_start ) ) then
    R_start = R
    Z_start = Z
  end if
end do
write(*,*) R_start, Z_start

n_bnd = 0
do i=1,element_list%n_elements  
  do iv = 1, n_vertex_max

    iv1 = iv
    iv2 = mod(iv,4) + 1

    inode1 = element_list%element(i)%vertex(iv1)
    inode2 = element_list%element(i)%vertex(iv2)

    if   (( ((node_list%node(inode1)%boundary .eq. 1) .or. (node_list%node(inode1)%boundary .eq. 3))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 1) .or. (node_list%node(inode2)%boundary .eq. 3)) ) &

!       .or. ( ((node_list%node(inode1)%boundary .eq. 2) .or. (node_list%node(inode1)%boundary .eq. 3))   &
!       .and.  ((node_list%node(inode2)%boundary .eq. 2) .or. (node_list%node(inode2)%boundary .eq. 3)) ) &

       .or. ( ((node_list%node(inode1)%boundary .eq. 4) .or. (node_list%node(inode1)%boundary .eq. 9))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 4) .or. (node_list%node(inode2)%boundary .eq. 9)) ) &

       .or. ( ((node_list%node(inode1)%boundary .eq. 4) .or. (node_list%node(inode1)%boundary .eq. 1))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 4) .or. (node_list%node(inode2)%boundary .eq. 1)) ) &

       .or. ( ((node_list%node(inode1)%boundary .eq. 5) .or. (node_list%node(inode1)%boundary .eq. 9))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 5) .or. (node_list%node(inode2)%boundary .eq. 9)) )) then 

      n_bnd = n_bnd + 1

    endif
  enddo
enddo

write(*,*) ' number of boundary points : ',n_bnd
write(*,*) ' number of toroidal planes : ',n_plane

nnos = n_plane * nsub * n_bnd
allocate(xyz(3,nnos),scalars(nnos,1:n_scalars),vectors(nnos,3,1:n_vectors))

nnoel = 4

if (periodic) then
  nel   = (n_plane)   * (nsub-1)*n_bnd
else
  nel   = (n_plane-1) * (nsub-1)*n_bnd
endif

allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0
n_points = nsub*n_bnd       ! number of points in one poloidal plane

allocate(avg6(n_points),avg7(n_points),avg8(n_points),avg11(n_points),Ravg(n_points),Zavg(n_points))
avg6 = 0.d0; avg7 = 0.d0; avg8 = 0.d0; avg11 = 0.d0; Ravg = 0.d0; Zavg = 0.d0

allocate(prf1(n_points),prf2(n_points),prf3(n_points),prf4(n_points),prf5(n_points),Rprf(n_points),Zprf(n_points))
prf2 = 0.d0; prf3 = 0.d0; prf4 = 0.d0; prf5 = 0.d0; Rprf = 0.d0; Zprf = 0.d0
allocate(prf6(n_points),prf7(n_points),prf8(n_points),prf11(n_points))
prf6 = 0.d0; prf7 = 0.d0; prf8 = 0.d0; prf11 = 0.d0

do i_tor=1, n_tor
  write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
enddo

do m=1,n_plane
  if (periodic) then
    phi = 2.d0 * PI * float(m-1)/float(n_plane)
  else
    phi = 2.d0 * PI * float(m-1)/float(n_plane-1) / float(n_period)
  endif
enddo

do m=1, n_plane

   inode_avg = 0

  if (periodic) then
    angle = 2.d0 * PI * float(m-1)/float(n_plane)
  else
    angle = 2.d0 * PI * float(m-1)/float(n_plane-1) / float(n_period)
  endif

  do i=1,element_list%n_elements

    do iv = 1, n_vertex_max

      iv1 = iv
      iv2 = mod(iv,4) + 1

      inode1 = element_list%element(i)%vertex(iv1)
      inode2 = element_list%element(i)%vertex(iv2)

      if (m .eq.1) then

        if ((node_list%node(inode1)%boundary .eq. 4) .and. (node_list%node(inode2)%boundary .eq. 5)) then
          write(*,*) ' PROBLEM mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 5) .and. (node_list%node(inode2)%boundary .eq. 4)) then
          write(*,*) ' PROBLEM2 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 1) .and. (node_list%node(inode2)%boundary .eq. 5)) then
          write(*,*) ' PROBLEM4 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 5) .and. (node_list%node(inode2)%boundary .eq. 1)) then
          write(*,*) ' PROBLEM6 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 2) .and. (node_list%node(inode2)%boundary .eq. 4)) then
          write(*,*) ' PROBLEM7 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 2) .and. (node_list%node(inode2)%boundary .eq. 5)) then
          write(*,*) ' PROBLEM8 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 4) .and. (node_list%node(inode2)%boundary .eq. 2)) then
          write(*,*) ' PROBLEM9 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif
        if ((node_list%node(inode1)%boundary .eq. 5) .and. (node_list%node(inode2)%boundary .eq. 2)) then
          write(*,*) ' PROBLEM10 mixed boundary',node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
        endif

      endif

      if   (( ((node_list%node(inode1)%boundary .eq. 1) .or. (node_list%node(inode1)%boundary .eq. 3))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 1) .or. (node_list%node(inode2)%boundary .eq. 3)) ) &
       .or. ( ((node_list%node(inode1)%boundary .eq. 4) .or. (node_list%node(inode1)%boundary .eq. 9))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 4) .or. (node_list%node(inode2)%boundary .eq. 9)) ) &
       .or. ( ((node_list%node(inode1)%boundary .eq. 4) .or. (node_list%node(inode1)%boundary .eq. 1))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 4) .or. (node_list%node(inode2)%boundary .eq. 1)) ) &
       .or. ( ((node_list%node(inode1)%boundary .eq. 5) .or. (node_list%node(inode1)%boundary .eq. 9))   &
       .and.  ((node_list%node(inode2)%boundary .eq. 5) .or. (node_list%node(inode2)%boundary .eq. 9)))) then 



        do j=1,nsub

         if   (  ((node_list%node(inode1)%boundary .eq. 1)  &
              .or.(node_list%node(inode1)%boundary .eq. 4)  &
              .or.(node_list%node(inode1)%boundary .eq. 9)  &
              .or.(node_list%node(inode1)%boundary .eq. 3)) &
           .and. ((node_list%node(inode2)%boundary .eq. 1)  &
              .or.(node_list%node(inode2)%boundary .eq. 4)  &
              .or.(node_list%node(inode2)%boundary .eq. 9)  &
              .or.(node_list%node(inode2)%boundary .eq. 3)) ) then
 
           sg = float(j-1)/float(nsub-1)

            if (iv1 .eq. 1) then
              tg     =  0.d0
              normal = -1.d0
            elseif (iv .eq. 3) then
              tg     = 1.d0
              normal = 1.d0
            else
              if ((m.eq.1) .and. (j.eq.1)) then
                write(*,*) ' problem 4/9 : ',i,iv1,iv2,inode1,inode2
                write(*,*) ' node1 R/Z   : ',node_list%node(inode1)%x(1,1,:),node_list%node(inode1)%boundary
                write(*,*) ' node2 R/Z   : ',node_list%node(inode2)%x(1,1,:),node_list%node(inode2)%boundary
              endif
            endif

         elseif (((node_list%node(inode1)%boundary .eq. 2) &
              .or.(node_list%node(inode1)%boundary .eq. 3)  &
              .or.(node_list%node(inode1)%boundary .eq. 5)  &
              .or.(node_list%node(inode1)%boundary .eq. 9)) &
           .and. ((node_list%node(inode2)%boundary .eq. 2)  &
              .or.(node_list%node(inode2)%boundary .eq. 3)  &
              .or.(node_list%node(inode2)%boundary .eq. 5)  &
              .or.(node_list%node(inode2)%boundary .eq. 9)) ) then

            tg = float(j-1)/float(nsub-1)

            if (iv1 .eq. 4) then
              sg     =  0.d0
              normal = -1.d0
            elseif (iv1 .eq.2) then
              sg     =  1.d0
              normal = +1.d0
            else
               if ((m .eq.1) .and. (j.eq.1)) write(*,*) ' problem 5/9 : ',i,iv1,iv2,inode1,inode2
            endif

          else
            if ((m .eq.1) .and. (j.eq.1)) then
              write(*,*) 'UNTREATED BOUNDARY : ',iv1,iv2,inode1,inode2,node_list%node(inode1)%boundary,node_list%node(inode2)%boundary
            endif
          endif

          call interp_RZ(node_list,element_list,i,sg,tg,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
          BigR = R
          xjac = R_s * Z_t - R_t * Z_s

          inode = inode+1

          inode_avg = inode_avg + 1

          xyz(1:3,inode) = (/ R * cos(angle), Z, R*sin(angle) /)

          Ravg(inode_avg) = R
          Zavg(inode_avg) = Z


          psi = 0.d0; psi_s = 0.d0; psi_t = 0.d0;
          rho = 0.d0; rho_s = 0.d0; rho_t = 0.d0;
          T   = 0.d0; T_s = 0.d0;   T_t = 0.d0; T_p = 0.d0
          Vpar = 0.d0

          do i_tor = 1,n_tor
            
            if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode

            call interp(node_list,element_list,i,1,i_tor,sg,tg,PS,PS_s,PS_t,PS_st,PS_ss,PS_tt)
            psi   = psi   + PS   * HZ(i_tor,m)
            psi_s = psi_s + PS_s * HZ(i_tor,m)
            psi_t = psi_t + PS_t * HZ(i_tor,m)

            call interp(node_list,element_list,i,5,i_tor,sg,tg,RH,RH_s,RH_t,RH_st,RH_ss,RH_tt)
            rho   = rho   + RH   * HZ(i_tor,m)
            rho_s = rho_s + RH_s * HZ(i_tor,m)
            rho_t = rho_t + RH_t * HZ(i_tor,m)

            call interp(node_list,element_list,i,6,i_tor,sg,tg,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
            T   = T   + TT   * HZ(i_tor,m)
            T_s = T_s + TT_s * HZ(i_tor,m)
            T_t = T_t + TT_t * HZ(i_tor,m)
            T_p = T_p + TT   * HZ_p(i_tor,m)

            call interp(node_list,element_list,i,7,i_tor,sg,tg,VP,VP_s,VP_t,VP_st,VP_ss,VP_tt)
            Vpar = Vpar + VP * HZ(i_tor,m)

          enddo
 
          psi_x = (   Z_t * psi_s - Z_s * psi_t ) / xjac
          psi_y = ( - R_t * psi_s + R_s * psi_t ) / xjac
          T_x   = (   Z_t * T_s   - Z_s * T_t )   / xjac
          T_y   = ( - R_t * T_s   + R_s * T_t )   / xjac

          BB2 = (F0**2 + (psi_x*psi_x+psi_y*psi_y)) / BigR**2
          
          psi_norm = get_psi_n(psi, Z)

          D_prof   = get_dperp (psi_norm)
          ZK_prof  = get_zkperp(psi_norm)

          ZKpar_T  = ZK_par * abs(max(T,T_min_ZKpar)/T_0)**2.5

          scalars(inode,1) = psi
          scalars(inode,2) = rho
          scalars(inode,3) = T
          scalars(inode,4) = Vpar * sqrt(BB2)

         if   (  ((node_list%node(inode1)%boundary .eq. 1)  &
              .or.(node_list%node(inode1)%boundary .eq. 4)  &
              .or.(node_list%node(inode1)%boundary .eq. 9)  &
              .or.(node_list%node(inode1)%boundary .eq. 3)) &
           .and. ((node_list%node(inode2)%boundary .eq. 1)  &
              .or.(node_list%node(inode2)%boundary .eq. 4)  &
              .or.(node_list%node(inode2)%boundary .eq. 9)  &
              .or.(node_list%node(inode2)%boundary .eq. 3)) ) then

            scalars(inode,5) = - (rho * Vpar * psi_s * normal)                  / R / sqrt(R_s**2 + Z_s**2)

            scalars(inode,6) = - (gamma_sheath -1.d0+gamma)*(rho * T * Vpar * psi_s * normal) / R / sqrt(R_s**2 + Z_s**2)

            if (abs(xjac) .gt. 1.d-7) then

            scalars(inode,7) = - ZKpar_T * (((T_x*psi_y - T_y*psi_x) + F0/BigR * T_p) / BigR / BB2 &
                             * (-psi_s*normal)) /BigR /sqrt(R_s**2+Z_s**2) 

            scalars(inode,8) = - ZK_prof * ( T_t   * (R_s**2 + Z_s**2) - T_s   * (R_s*R_t+Z_s*Z_t)) / Xjac &
                             / sqrt(R_s**2+Z_s**2) * normal

            scalars(inode,9) = - D_prof  * ( rho_t * (R_s**2 + Z_s**2) - rho_s * (R_s*R_t+Z_s*Z_t)) / Xjac &
                             / sqrt(R_s**2+Z_s**2) * normal
            endif

            scalars(inode,10) = - psi_s / BigR /sqrt(R_s**2+Z_s**2) * normal
            scalars(inode,11) = gamma_sheath * rho * T * Vpar * sqrt(BB2) 
            scalars(inode,12) = normal
            scalars(inode,13) = - rho * (Vpar * sqrt(BB2))**2 * Vpar * psi_s * normal / R / sqrt(R_s**2 + Z_s**2)

            vectors(inode,:,1) = (/ + psi_y /BigR * cos(angle), - psi_x /BigR, + psi_y /BigR * sin(angle) /) 
            vectors(inode,:,2) = (/ + vpar * psi_y /BigR* cos(angle), - vpar * psi_x /BigR, + vpar * psi_y /BigR * sin(angle) /) 
            vectors(inode,:,3) = (/ - Z_s * cos(angle), + R_s, -Z_s * sin(angle)  /) / sqrt(R_s**2+Z_s**2) * normal

            avg6(inode_avg) = avg6(inode_avg) + scalars(inode,6)
            avg7(inode_avg) = avg7(inode_avg) + scalars(inode,7)
            avg8(inode_avg) = avg8(inode_avg) + scalars(inode,8)
            avg11(inode_avg) = avg11(inode_avg) + scalars(inode,11)

         elseif (((node_list%node(inode1)%boundary .eq. 2) &
              .or.(node_list%node(inode1)%boundary .eq. 3)  &
              .or.(node_list%node(inode1)%boundary .eq. 5)  &
              .or.(node_list%node(inode1)%boundary .eq. 9)) &
           .and. ((node_list%node(inode2)%boundary .eq. 2)  &
              .or.(node_list%node(inode2)%boundary .eq. 3)  &
              .or.(node_list%node(inode2)%boundary .eq. 5)  &
              .or.(node_list%node(inode2)%boundary .eq. 9)) ) then

            scalars(inode,5) = (rho * Vpar * psi_t * normal)                    / R / sqrt(R_t**2 + Z_t**2)

            scalars(inode,6) = gamma_sheath * (rho * T * Vpar * psi_t * normal) / R / sqrt(R_t**2 + Z_t**2)
 
            if (abs(xjac) .gt. 1.d-7) then

            scalars(inode,7) = - ZKpar_T * (((T_x*psi_y - T_y*psi_x) + F0/BigR * T_p) / BigR / BB2 &
                             * (psi_t*normal)) /BigR /sqrt(R_t**2+Z_t**2) 

            scalars(inode,8) = - ZK_prof * ( T_s   * (R_t**2 + Z_t**2) - T_t   * (R_s*R_t+Z_s*Z_t)) / Xjac &
                             / sqrt(R_t**2+Z_t**2) * normal

            scalars(inode,9) = - D_prof  * ( rho_s * (R_t**2 + Z_t**2) - rho_t * (R_s*R_t+Z_s*Z_t)) / Xjac &
                             / sqrt(R_t**2+Z_t**2) * normal

            endif

            scalars(inode,10) = psi_t / BigR /sqrt(R_t**2+Z_t**2) * normal
            scalars(inode,11) = gamma_sheath * rho * T * Vpar * sqrt(BB2) 
            scalars(inode,12) = normal
            scalars(inode,13) = -rho * (Vpar * sqrt(BB2))**2 * Vpar * psi_t * normal / R / sqrt(R_t**2 + Z_t**2)

            vectors(inode,:,1) = (/ + psi_y /BigR * cos(angle), - psi_x /BigR, + psi_y /BigR * sin(angle) /) 
            vectors(inode,:,2) = (/ + vpar * psi_y /BigR* cos(angle), - vpar * psi_x /BigR, + vpar * psi_y /BigR * sin(angle)  /) 
            vectors(inode,:,3) = (/ + Z_t * cos(angle), - R_t, Z_t * sin(angle)  /) / sqrt(R_t**2+Z_t**2) * normal

            avg6(inode_avg) = avg6(inode_avg) + scalars(inode,6)
            avg7(inode_avg) = avg7(inode_avg) + scalars(inode,7)
            avg8(inode_avg) = avg8(inode_avg) + scalars(inode,8)
            avg11(inode_avg) = avg11(inode_avg) + scalars(inode,11)

          else

            write(*,*) 'warning untreated boundary : ',inode1,inode2,node_list%node(inode1)%boundary,node_list%node(inode2)%boundary

          endif

          if (m.eq.1) then
            Rprf(inode_avg)  = R
            Zprf(inode_avg)  = Z
	    prf1(inode_avg)  = scalars(inode,1)
            prf7(inode_avg)  = scalars(inode,7) / MU_zero / t_norm * 1.5
            prf6(inode_avg)  = scalars(inode,6) / MU_zero / t_norm * 1.5
            prf2(inode_avg)  = scalars(inode,2) * central_density
            prf3(inode_avg)  = scalars(inode,3) / MU_zero / (central_density * 1d20) / 1.602d-19 /2.
            prf4(inode_avg)  = scalars(inode,4) / t_norm
            prf5(inode_avg)  = scalars(inode,5) * central_density / t_norm
            prf11(inode_avg) = scalars(inode,11) / MU_zero / t_norm
          endif

        enddo

        if (m .lt. n_plane) then

          do j=1,nsub-1

            ielm        = ielm + 1
            ien(1,ielm) = inode - nsub + (j-1)      ! 0 based indices for VTK
            ien(2,ielm) = inode - nsub + (j  ) 
            ien(3,ielm) = ien(2,ielm) + n_points
            ien(4,ielm) = ien(1,ielm) + n_points

          enddo

        endif

        if ( (periodic) .and. (m .eq. n_plane)) then

          do j=1,nsub-1

            ielm        = ielm+1
            ien(1,ielm) = inode - nsub + (j-1)       ! 0 based indices for VTK
            ien(2,ielm) = inode - nsub + (j  ) 

            ien(3,ielm) = ien(2,ielm) - n_points * (n_plane-1)
            ien(4,ielm) = ien(1,ielm) - n_points * (n_plane-1)

          enddo

        endif
	
      endif

    enddo
  enddo
enddo

close(22)


!scalar_names = (/ 'flux    ','density ','T       ','Vpar    ','nV.n    ','nTV.n   ','KparT.n ', &
!                  'Kperp.T ','Dperp.n ','B.n     '/)

scalars(:,2) = scalars(:,2) * central_density
scalars(:,3) = scalars(:,3) / MU_zero / (central_density * 1d20) / 1.602d-19 /2. !(assumes Te=Ti=T/2)
scalars(:,4) = scalars(:,4) / t_norm
scalars(:,5) = scalars(:,5) * central_density / t_norm
scalars(:,6) = scalars(:,6) / MU_zero / t_norm * 1.5
scalars(:,7) = scalars(:,7) / MU_zero / t_norm * 1.5
scalars(:,8) = scalars(:,8) / MU_zero / t_norm * 1.5
scalars(:,9) = scalars(:,9) * central_density / t_norm
scalars(:,11) = scalars(:,11) / MU_zero / t_norm
scalars(:,13) = scalars(:,13) / MU_zero / t_norm * 0.5

avg6(:) = avg6(:) / MU_zero / t_norm / float(n_plane) * 1.5
avg7(:) = avg7(:) / MU_zero / t_norm / float(n_plane) * 1.5
avg8(:) = avg8(:) / MU_zero / t_norm / float(n_plane) * 1.5
avg11(:)= avg11(:) / MU_zero / t_norm / float(n_plane)


open(22,file='target_profile')
write(22,'(A132)') "  Length          R               Z               angle           KparT_normal    gam_nVT_normal  density         T              Vpar            nv_normal       gam_nvT"

open(23,file='average_target_profile')
write(23,'(A132)') '      time         step            Length         R               Z              nTV.n           KparT.n        Kperp.T        nvT_gam'

allocate(tobedone(n_points))

tobedone = .true.
distance = 0.d0

do i=1, n_points

  index = minval(minloc((Rprf-R_start)**2+(Zprf-Z_start)**2,tobedone))

  distance = distance + sqrt((Rprf(index)-R_start)**2+(Zprf(index)-Z_start)**2)

  R_start = Rprf(index)
  Z_start = Zprf(index)
  tobedone(index) = .false.

  write(22,'(12e16.8)') distance, Rprf(index),Zprf(index),angle,prf7(index),prf6(index), &
                        prf2(index), prf3(index),  prf4(index), prf5(index), prf11(index) 

  write(23,'(12e16.8)') xtime(index_start),xtime(index_start)-xtime(index_start-10), distance, &
                        Ravg(index),Zavg(index),avg6(index),avg7(index),avg8(index),avg11(index),normal
enddo

close(23)
close(22)

!--------------------------------------------------- write the binary VTK file
etype = 9  ! for vtk_quad

lf = char(10) ! line feed character

open(unit=ivtk,file='jorek_tmp.vtk',access='stream',form='unformatted',convert='BIG_ENDIAN',status='replace')

buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)

! POINTS SECTION
write(str1(1:8),'(i8)') nnos
buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
write(ivtk) ((xyz(i,j),i=1,3),j=1,nnos)

! CELLS SECTION
write(str1(1:8),'(i8)') nel            ! number of elements (cells)
write(str2(1:8),'(i8)') nel*(1+nnoel)  ! size of the following element list (nel*(nnoel+1))
buffer = lf//lf//'CELLS '//str1//' '//str2//lf                                        ; write(ivtk) trim(buffer)
write(ivtk) (nnoel,(ien(i,j),i=1,nnoel),j=1,nel)

! CELL_TYPES SECTION
write(str1(1:8),'(i8)') nel   ! number of elements (cells)
buffer = lf//lf//'CELL_TYPES'//str1//lf                                               ; write(ivtk) trim(buffer)
write(ivtk) (etype,i=1,nel)

! POINT_DATA SECTION
write(str1(1:8),'(i8)') nnos
buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)

do i_var =1, n_scalars
  buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
  write(ivtk) (scalars(i,i_var),i=1,nnos)
enddo

do i_var =1, n_vectors
  buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf                      ; write(ivtk) trim(buffer)
  write(ivtk) ((vectors(j,i,i_var),i=1,3),j=1,nnos)
enddo

close(ivtk)

end


