!> Program to convert a JOREK2 restart file into volume/surface averaged Vorticity Terms and VTK file
!---------------------------------------------------------------------------------------------------
!
! To compile:
!   * make jorek2vtk_GaussVortTerms
! 
! The executable "jorek2vtk_GaussVortTerms" is created.
!
! INPUT:
!
! The executable needs the following files:
!   * jorek_restart.rst|h5 (main jorek calculation data)
!   * the INPUT file use for the main jorek calculation
!   * (OPTIONAL) vtk_GaussVortTerms.nml (in this file: if n_plane_local = n_plane -> volume averaged 
!      calculation, if n_plane_local = 1 -> surface averaged calculation)
!
! EXECUTE:
!
! To execute for a single time step (for the current jorek_restart.rst|h5):
!   * ./jorek2vtk_GaussVortTerms < INPUT_file
!
! To run for several time steps use the bash script: run_vtk_GaussTerms.sh
! Example: if you want to execute from the iteration 1000 to 3000 by steps of 20, 
!          first check the file run_vtk_GaussTerms.sh and change the INPUT file name on line 42,
!          then run:
!            * ./run_vtk_GaussTerms.sh 1000 3000 20
! 
! OUTPUT:
!
! Files created after the run (FOR A SINGLE TIME STEP):
!   * Average_Terms_Voriticity_Gauss.dat (Volume/surface averaged terms of the vorticity equation)
!   * jorek_tmp.vtk                      (vtk file, to be read with paraview or other)
!   * Diff_Vort_Terms_in_Element.dat     (Value of the different vorticity terms for all the elements
!                                         calculated in weak form with integration by parts and 
!                                         noted with "_2" at the end if calculated WITHOUT integration
!                                         by parts) 
!
! Files created after the run (FOR MULTIPLE TIME STEPS, using the script run_vtk_GaussTerms.sh):
!
!   * Average_Terms_Voriticity_Gauss.dat    (Volume/surface averaged terms of the vorticity equation)
!   * Average_Terms_Voriticity_Gauss_IT.dat (Volume/surface averaged terms of the vorticity equation
!                                            at the different times/ITERATIONS)
!   * jorekGaussXXXXX.vtk                   (vtk files at ITERATION, to be read with paraview or other)
!   * Diff_Vort_Terms_in_Element.dat        (Value of the different vorticity terms for all the elements
!                                            calculated in weak form with integration by parts and 
!                                            noted with "_2" at the end if calculated WITHOUT integration
!                                            by parts, JUST FOR THE LAST ITERATION CONSIDERED)
!
!---------------------------------------------------------------------------------------------------

! --- PROGRAM
program jorek2vtk_GaussVortTerms
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use corr_neg
use elements_nodes_neighbours
use mod_import_restart
use mod_neighbours
use equil_info, only : get_psi_n, ES
use mod_sources


implicit none

!type (type_node_list)    :: node_list
!type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

integer    :: i, j, k, in, ms, mt, mp, im, iv, inode, ife, n_elements
integer    :: k_tor, ierr, my_id, i_elm_axis, i_elm_xpoint, ifail, n_plane_local

real*8     :: psi_norm
real*8     :: wst, xjac, xjac_x, xjac_y, error, Btot, BigR, BigZ

real*8     :: current_source(n_gauss,n_gauss), particle_source(n_gauss,n_gauss)                    &
            , heat_source(n_gauss,n_gauss)

real*8, dimension(n_plane,n_var,n_gauss,n_gauss)  :: eq_g, eq_s, eq_t
real*8, dimension(n_plane,n_var,n_gauss,n_gauss)  :: eq_p
real*8, dimension(n_plane,n_var,n_gauss,n_gauss)  :: eq_ss, eq_st, eq_tt   
real*8, dimension(n_plane,n_var,n_gauss,n_gauss)  :: delta_g, delta_s, delta_t, delta_ss, delta_tt &
                                                   , delta_st

real*8, dimension(n_plane,n_gauss,n_gauss)        :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt  &
                                                   , v_xx, v_yy, v_xy

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8     :: delta_u, delta_u_s, delta_u_t, delta_u_ss, delta_u_tt, delta_u_st

real*8     :: ZK_prof, D_prof, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y, ZKpar_T  &
            , dZKpar_dT, eps_cyl, delta_u_xx, delta_u_yy

real*8     :: ps0, ps0_x, ps0_y, ps0_p,ps0_s,ps0_t, ps0_ss, ps0_st, ps0_tt, ps0_xx, ps0_xy, ps0_yy
real*8     :: zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_st, u0_tt, u0_xx, u0_xy, u0_yy 
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt     &
            , vpar0_xx, vpar0_xy, vpar0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy, r0_hat   &
            , r0_x_hat, r0_y_hat
real*8     :: T0, T0_x, T0_y, T0_p, T0_s, T0_t, T0_ss, T0_st, T0_tt, T0_xx, T0_xy, T0_yy
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_yy, psi_xy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t, zj_ss, zj_st, zj_tt
real*8     :: vpar, vpar_x, vpar_y, vpar_s, vpar_t, vpar_p, vpar_ss, vpar_st, vpar_tt, vpar_xx     &
            , vpar_xy, vpar_yy
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy
real*8     :: rho_hat, rho_x_hat, rho_y_hat
real*8     :: T, T_x, T_y, T_s, T_t, T_p, T_ss, T_st, T_tt, T_xx, T_xy, T_yy
real*8     :: P0, P0_x, P0_y, P0_s, P0_t, P0_ss, P0_st, P0_tt, P0_p, P0_xx, P0_xy, P0_yy
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T, eta_num_T
real*8     :: psi_abs, Btheta

character(len=512), parameter :: MODE_write = "(3x,I3,': ',A,'(',A,'*phi)')"
character(len=10)             :: mode_num_TermSum

! --- Diag Terms vorticity equation
logical               :: exist_file
integer, parameter    :: file_handle = 99999
integer               :: errorFile, n_tor_Terms, n_period_Terms, n_plane_Terms
character(len=1024)   :: filename
real*8                :: rot_NL_Ve, gd_rho_NL_Ve, Poiss_Psi_j, der_j_tor                   &
                       , rot_grad_P, rot_Vstar, visco_lap_W, Source_term, rot_delta_u      &
                       , visco_num_W
real*8                :: Sum_rot_NL_Ve, Sum_gd_rho_NL_Ve, Sum_Poiss_Psi_j, Sum_der_j_tor           &
                       , Sum_rot_grad_P, Sum_rot_Vstar, Sum_visco_lap_W, Sum_Source_term           &
                       , Sum_vv2, Sum_rot_delta_u, Sum_visco_num_W, Sum_Surface, Sum_rot_grad_P_2  &
                       , Sum_xjac, Sum_sum_1, Sum_Sum_2
integer               :: count_Ve, count_Ve2

real*8                :: RHS_local(44,1,1)
real*8                :: Poiss_Psi_j_2, rot_grad_P_2, rot_NL_Ve_2, gd_rho_NL_Ve_2, visco_lap_W_2   &
                       , rot_delta_u_2, Source_term_2, w0_2xx, w0_2yy, visco_num_W_2, rot_Vstar_2

! --- For VTK files
real*4 , allocatable        :: xyz(:,:), scalars(:,:)
integer, allocatable        :: ien(:,:), ien2(:,:)
integer                     :: n_scalars, n_terms, etype, i_var, iside_i, iside_j                  &
                             , n_element_less_private
character*12, allocatable   :: scalar_names(:)
integer, parameter          :: ivtk = 22 ! an arbitrary unit number for the VTK output file
character                   :: buffer*80, lf*1, str1*12, str2*12

namelist /vtk_params/          n_plane_local 



! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Initialize mode and mode_type arrays
call det_modes()

! --- Read rst file
call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

! --- Define the basis functions at the Gaussian points
call initialise_basis

! --- Write modes used
do k_tor = 2, n_tor
  write(mode_num_TermSum,'(I4)') mode(k_tor)
  write(*,MODE_write) k_tor, mode_type(k_tor), trim(adjustl(mode_num_TermSum))
end do

! --- Initialization for VTK files
n_terms   = 22
n_scalars = n_var + n_terms

allocate(scalar_names(n_scalars))

scalar_names(1:n_var) = variable_names((/(i, i=1,n_var)/))

scalar_names(n_var+1:n_var+n_terms) = (/                                   &
    'Er          ' , 'Vtheta      ' , 'Vsound      ' , 'Btot_T      '      &
  , 'psi_norm    ' , 'Vt_ExB      ' , 'Vt_Star     ' , 'Vt_V0       '      &
  , 'rotDelta_Vt ' , 'Poiss_Psi_j ' , 'rot_grad_P  ' , 'Dj_D_tor    '      &
  , 'rot_NL_Ve   ' , 'Grd_rho_NL  ' , 'Visco_W     ' , 'Visco_NUM_W '      &
  , 'Src_rho_V   ' , 'rt_Vstar_NL ' , 'VV2         ' , 'Equilibrium '      &
  , 'Equil+Visco ' , 'Ptc_Source  '  /)
  
allocate( xyz(3,element_list%n_elements),scalars(element_list%n_elements,1:n_scalars) )
allocate( ien (4,(element_list%n_elements)) )
allocate( ien2(4,(element_list%n_elements)) )
xyz     = 0.
scalars = 0.
ien     = 0
ien2    = 0

! --- Volume average of the vorticity terms equation; initialization:
Sum_rot_NL_Ve    = 0.
Sum_gd_rho_NL_Ve = 0.
Sum_Poiss_Psi_j  = 0.
Sum_der_j_tor    = 0.
Sum_rot_grad_P   = 0.
Sum_rot_grad_P_2 = 0.
Sum_visco_lap_W  = 0.
Sum_visco_num_W  = 0.
Sum_rot_Vstar    = 0.
Sum_Source_term  = 0.
Sum_vv2          = 0.
Sum_rot_delta_u  = 0.
Sum_Surface      = 0.
Sum_xjac         = 0.
Sum_sum_1        = 0.
Sum_sum_2        = 0.

rot_NL_Ve        = 0.
gd_rho_NL_Ve     = 0.
Poiss_Psi_j      = 0.
der_j_tor        = 0.
rot_grad_P       = 0.
visco_lap_W      = 0.
rot_Vstar        = 0.
Source_term      = 0.
vv2              = 0.
rot_delta_u      = 0. 

count_Ve  = 0
count_Ve2 = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! --- Maximum value for n_plane_local must be n_plane and minimum is 1 !!!
n_plane_local = n_plane ! --- Value by default

! --- Read parameters from namelist file 'vtk_GaussVortTerms.nml' if it exists
open(42, file='vtk_GaussVortTerms.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from vtk_GaussVortTerms.nml namelist.'
  read(42,vtk_params)
  close(42)
  
  if ( n_plane_local .GT. n_plane ) then
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print*, ' ERROR  n_plane_local  MUST be LESS OR EQUAL than  n_plane !!! '
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    STOP
  endif
else
  print*, ' '
  print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*, 'NO vtk_GaussVortTerms.nml file using : n_plane_local = n_plane , by default!'
  print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*, ' '
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*) '                         '
write(*,*) 'Parameters:'
write(*,*) '-----------'
write(*,*) 'n_tor                  = ', n_tor
write(*,*) 'n_period               = ', n_period
write(*,*) 'n_plane                = ', n_plane
write(*,*) 'n_plane_local          = ', n_plane_local
write(*,*) 'n_order                = ', n_order
write(*,*) 'n_vertex_max           = ', n_vertex_max
write(*,*) 'n_gauss                = ', n_gauss
write(*,*) 'n_var                  = ', n_var
write(*,*) 'n_scalars              = ', n_scalars
write(*,*) 'n_tht                  = ', n_tht
write(*,*) 'n_flux                 = ', n_flux
write(*,*) 'n_open                 = ', n_open
write(*,*) 'n_leg                  = ', n_leg
write(*,*) 'n_private              = ', n_private
write(*,*) 'tstep_n                = ', tstep_n
write(*,*) 'wgauss                 = ', wgauss
write(*,*) 'visco                  = ', visco
write(*,*) 'visco_num              = ', visco_num
write(*,*) 'xpoint                 = ', xpoint
write(*,*) 'visco_T_dependent      = ', visco_T_dependent
write(*,*) 'psi_axis               = ', ES%psi_axis
write(*,*) 'psi_bnd                = ', ES%psi_bnd
write(*,*) 'psi_xpoint             = ', ES%psi_xpoint
write(*,*) 'Z_xpoint               = ', ES%Z_xpoint
write(*,*) 'F0                     = ', F0
write(*,*) 'R_geo,Z_geo            = ', R_geo, Z_geo
write(*,*) 'HZ                     = ', HZ
write(*,*) '                         '
write(*,*) 'particlesource         = ', particlesource
write(*,*) 'particlesource_psin    = ', particlesource_psin
write(*,*) 'particlesource_sig     = ', particlesource_sig
write(*,*) '                         '


! --- Write the different terms of the vorticity equation for the current time step at the 
! --- center of the element (mean R, Z value) point in the poloidal plane
open(unit=999992, file='Diff_Vort_Terms_in_Element.dat', status='unknown',action='write'  &
   , iostat=errorFile)

if (errorFile /= 0) then
  write(*,*) 'ERROR : Opening file Diff_Vort_Terms_in_Element"' &
            , trim(filename), '" failed.'
end if

write(999992, '(1A7,30A16)')  , ' n1 ife'           &
                              , ' n2 BigR        '  &
                              , ' n3 BigZ        '  & 
                              , ' n4 ps0        '   &
                              , ' n5 rotDelta_u  '  &
                              , ' n6 rotDelta_u_2'  &
                              , ' n7 PoissPsi_j  '  &
                              , ' n8 PoissPsi_j_2'  &
                              , ' n9 rot_gd_P    '  &
                              , ' n10 rot_gd_P_2 '  &
                              , ' n11 der_j_tor  '  &
                              , ' n12 rot_NL_Ve  '  &
                              , ' n13 rot_NL_Ve_2'  &
                              , ' n14 rho_NL_Ve  '  &
                              , ' n15 rho_NL_Ve_2'  &
                              , ' n16 visco_W    '  &
                              , ' n17 visco_W_2  '  &
                              , ' n18 viscoNumW  '  &
                              , ' n19 viscoNumW_2'  &
                              , ' n20 SrceTerm   '  &
                              , ' n21 SrceTerm_2 '  &
                              , ' n22 v_s        '  &
                              , ' n23 v_t        '  &
                              , ' n24 xjac       '  &
                              , ' n25 rot_Vstar  '  &
                              , ' n26 rot_Vstar_2'  &
                              , ' n27 r0 * vv2   '  &
                              , ' n28 SumEqui    '  &
                              , ' n29 SumEquiVis '  &
                              , ' n30 Mask(1or 0)'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   BEGIN LOOP ELEMENTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ife = 1, element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

 !--- Value of (x,y) and derivatives on Gaussian points, initialization
 x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0
 y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0
 eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0; eq_p = 0.d0
 
 ! --- Test function initialization 
 v    = 0.d0; v_x  = 0.d0; v_y  = 0.d0; v_s  = 0.d0; v_t = 0.d0; v_p = 0.d0; v_ss = 0.d0
 v_st = 0.d0; v_tt = 0.d0; v_xx = 0.d0; v_yy = 0.d0; v_xy = 0.d0 
 
 !print*, 'v             = ', v
 !print*, 'v_s           = ', v_s
 !print*, 'v_t           = ', v_t

 ! --- Delta initialization 
 delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0; delta_ss = 0.d0; delta_tt = 0.d0; delta_st = 0.d0

 do i=1,n_vertex_max
  do j=1,n_order+1
   do ms=1,n_gauss
    do mt=1,n_gauss

      x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
      x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
      x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

      x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
      x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
      x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

      y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
      y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
      y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

      y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
      y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)
      y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)
      
      do mp=1,n_plane_local
       do k=1,n_var
        do in=1,n_tor

          eq_g(mp,k,ms,mt)    = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j)   &
                              * H(i,j,ms,mt)  * HZ(in,mp)
          eq_s(mp,k,ms,mt)    = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j)   &
                              * H_s(i,j,ms,mt)* HZ(in,mp)
          eq_t(mp,k,ms,mt)    = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j)   &
                              * H_t(i,j,ms,mt)* HZ(in,mp)
          eq_p(mp,k,ms,mt)    = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j)   &
                              * H(i,j,ms,mt)  * HZ_p(in,mp)

          eq_ss(mp,k,ms,mt)   = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j)  &
                              * H_ss(i,j,ms,mt) * HZ(in,mp)
          eq_st(mp,k,ms,mt)   = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j)  &
                              * H_st(i,j,ms,mt) * HZ(in,mp)
          eq_tt(mp,k,ms,mt)   = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j)  &
                              * H_tt(i,j,ms,mt) * HZ(in,mp)

          delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k)  &
                              * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
          delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k)  &
                              * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
          delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k)  &
                              * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)
                              
          delta_ss(mp,k,ms,mt) = delta_ss(mp,k,ms,mt) + nodes(i)%deltas(in,j,k)  &
                               * element%size(i,j) * H_ss(i,j,ms,mt) * HZ(in,mp)
          delta_st(mp,k,ms,mt) = delta_st(mp,k,ms,mt) + nodes(i)%deltas(in,j,k)  &
                               * element%size(i,j) * H_st(i,j,ms,mt) * HZ(in,mp)
          delta_tt(mp,k,ms,mt) = delta_tt(mp,k,ms,mt) + nodes(i)%deltas(in,j,k)  &
                               * element%size(i,j) * H_tt(i,j,ms,mt) * HZ(in,mp)

        enddo
       enddo
      enddo
       
    enddo
   enddo
  enddo
 enddo
 
 do i=1,n_vertex_max
  do j=4,4  ! --- Just the 4th type test function is used because is the only to be zero at boundaries
   do ms=1,n_gauss
    do mt=1,n_gauss
 
      xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
                    
      if (abs(xjac) .LT. 1E-15) print*, 'xjac = ', xjac

      do mp=1,n_plane_local
       do im=1,1            ! --- Axisymmetric test function
       
         v(mp,ms,mt)   = v(mp,ms,mt)   + H(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)

         v_x(mp,ms,mt) = v_x(mp,ms,mt) + (  y_t(ms,mt) * h_s(i,j,ms,mt)                &
                                          - y_s(ms,mt) * h_t(i,j,ms,mt) )              &
                       * element%size(i,j) / xjac * HZ(im,mp)

         v_y(mp,ms,mt) = v_y(mp,ms,mt) + (- x_t(ms,mt) * h_s(i,j,ms,mt)                &
                                          + x_s(ms,mt) * h_t(i,j,ms,mt) )              &
                       * element%size(i,j) / xjac * HZ(im,mp)

         v_s(mp,ms,mt) = v_s(mp,ms,mt) + h_s(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
         v_t(mp,ms,mt) = v_t(mp,ms,mt) + h_t(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
         v_p(mp,ms,mt) = v_p(mp,ms,mt) + H(i,j,ms,mt)   * element%size(i,j) * HZ_p(im,mp)

         v_ss(mp,ms,mt) = v_ss(mp,ms,mt) + h_ss(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
         v_tt(mp,ms,mt) = v_tt(mp,ms,mt) + h_tt(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
         v_st(mp,ms,mt) = v_st(mp,ms,mt) + h_st(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
         
       enddo
      enddo
      
    enddo
   enddo
  enddo
 enddo
 
 ! --- Test function derivatives in global coordinates calculation
 do ms=1,n_gauss
  do mt=1,n_gauss
      
    xjac    = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
    
    if (abs(xjac) .LT. 1E-15) print*, 'xjac = ', xjac
      
    xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt)         &
            - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)                                 &
            + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))            &
            + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

    xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt)         &
            - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)                                 &
            + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))            &
            + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac
              
    do mp=1,n_plane_local

      v_xx(mp,ms,mt) = ( v_ss(mp,ms,mt) * y_t(ms,mt)**2                                            &
                     -   2.d0*v_st(mp,ms,mt) * y_s(ms,mt) * y_t(ms,mt)                             &
                     +   v_tt(mp,ms,mt) * y_s(ms,mt)**2                                            &
                     +   v_s(mp,ms,mt)  * (y_st(ms,mt) * y_t(ms,mt) - y_tt(ms,mt) * y_s(ms,mt) )   &
                     +   v_t(mp,ms,mt)  * (y_st(ms,mt) * y_s(ms,mt) - y_ss(ms,mt) * y_t(ms,mt) ) ) &
                         / xjac**2                                                                 &
                     -   xjac_x * (   v_s(mp,ms,mt) * y_t(ms,mt) - v_t(mp,ms,mt) * y_s(ms,mt) )    &
                         / xjac**2

      v_yy(mp,ms,mt) = ( v_ss(mp,ms,mt) * x_t(ms,mt)**2                                            &
                     -   2.d0*v_st(mp,ms,mt) * x_s(ms,mt) * x_t(ms,mt)                             &
                     +   v_tt(mp,ms,mt) * x_s(ms,mt)**2                                            &
                     +   v_s(mp,ms,mt)  * (x_st(ms,mt) * x_t(ms,mt) - x_tt(ms,mt) * x_s(ms,mt) )   &
                     +   v_t(mp,ms,mt)  * (x_st(ms,mt) * x_s(ms,mt) - x_ss(ms,mt) * x_t(ms,mt) ) ) &
                         / xjac**2                                                                 &
                     -   xjac_y * ( - v_s(mp,ms,mt) * x_t(ms,mt) + v_t(mp,ms,mt) * x_s(ms,mt) )    &
                         / xjac**2

      v_xy(mp,ms,mt) = ( - v_ss(mp,ms,mt) * y_t(ms,mt) * x_t(ms,mt)                                &
                     -     v_tt(mp,ms,mt) * x_s(ms,mt) * y_s(ms,mt)                                &
                     +     v_st(mp,ms,mt) * (y_s(ms,mt)  * x_t(ms,mt) + y_t(ms,mt)  * x_s(ms,mt))  &
                     -     v_s(mp,ms,mt)  * (x_st(ms,mt) * y_t(ms,mt) - x_tt(ms,mt) * y_s(ms,mt))  &
                     -     v_t(mp,ms,mt)  * (x_st(ms,mt) * y_s(ms,mt) - x_ss(ms,mt) * y_t(ms,mt)) )&
                           / xjac**2                                                               &
                     -     xjac_x * (- v_s(mp,ms,mt) * x_t(ms,mt) + v_t(mp,ms,mt) * x_s(ms,mt) )   &
                          / xjac**2
    enddo
  enddo
 enddo
 
 ! --- Initialize table RHS_local to zero for each element
 RHS_local = 0.
 
 ! --- Sum over the Gaussian integration points
 do ms=1, n_gauss
  do mt=1, n_gauss

    wst = wgauss(ms)*wgauss(mt)

    xjac    = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
    

    xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt)         &
            -  2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)                                &
            +  y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))           &
            +  x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

    xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt)         &
            -  2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)                                &
            +  x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))           &
            +  y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

    BigR    = x_g(ms,mt)
    BigZ    = y_g(ms,mt)
    BigR_x  = 1.d0

    eps_cyl = 1.d0          ! for cylinder geometry : epscyl = eps

#ifdef WITH_TiTe
    write(*,*) 'ERROR: Diagnostic has not been adapted to models with separate Te and Ti.'
    stop
#else
    call sources(xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(1,1,ms,mt), ES%psi_axis, ES%psi_bnd  &
               , particle_source(ms,mt), heat_source(ms,mt))
#endif

    do mp = 1, n_plane_local

      ! --- Variables and variables derivatives definition
      ps0    = eq_g(mp,1,ms,mt)
      ps0_x  = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
      ps0_y  = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
      ps0_p  = eq_p(mp,1,ms,mt)
      ps0_s  = eq_s(mp,1,ms,mt)
      ps0_t  = eq_t(mp,1,ms,mt)
      ps0_ss = eq_ss(mp,1,ms,mt)
      ps0_tt = eq_tt(mp,1,ms,mt)
      ps0_st = eq_st(mp,1,ms,mt)

      u0    = eq_g(mp,2,ms,mt)
      u0_x  = (   y_t(ms,mt) * eq_s(mp,2,ms,mt) - y_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
      u0_y  = ( - x_t(ms,mt) * eq_s(mp,2,ms,mt) + x_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
      u0_p  = eq_p(mp,2,ms,mt)
      u0_s  = eq_s(mp,2,ms,mt)
      u0_t  = eq_t(mp,2,ms,mt)
      u0_ss = eq_ss(mp,2,ms,mt)
      u0_tt = eq_tt(mp,2,ms,mt)
      u0_st = eq_st(mp,2,ms,mt)

      vv2   = BigR**2 *  ( u0_x * u0_x + u0_y *u0_y  )

      zj0   = eq_g(mp,3,ms,mt)
      zj0_x = (   y_t(ms,mt) * eq_s(mp,3,ms,mt) - y_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
      zj0_y = ( - x_t(ms,mt) * eq_s(mp,3,ms,mt) + x_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
      zj0_p = eq_p(mp,3,ms,mt)
      zj0_s = eq_s(mp,3,ms,mt)
      zj0_t = eq_t(mp,3,ms,mt)

      w0    = eq_g(mp,4,ms,mt)
      w0_x  = (   y_t(ms,mt) * eq_s(mp,4,ms,mt) - y_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
      w0_y  = ( - x_t(ms,mt) * eq_s(mp,4,ms,mt) + x_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
      w0_p  = eq_p(mp,4,ms,mt)
      w0_s  = eq_s(mp,4,ms,mt)
      w0_t  = eq_t(mp,4,ms,mt)
      w0_ss = eq_ss(mp,4,ms,mt)
      w0_tt = eq_tt(mp,4,ms,mt)
      w0_st = eq_st(mp,4,ms,mt)

      r0    = eq_g(mp,5,ms,mt)
      r0_x  = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
      r0_y  = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
      r0_p  = eq_p(mp,5,ms,mt)
      r0_s  = eq_s(mp,5,ms,mt)
      r0_t  = eq_t(mp,5,ms,mt)
      r0_ss = eq_ss(mp,5,ms,mt)
      r0_st = eq_st(mp,5,ms,mt)
      r0_tt = eq_tt(mp,5,ms,mt)

      r0_hat   = BigR**2 * abs(r0)
      r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
      r0_y_hat = BigR**2 * r0_y

      T0    = eq_g(mp,6,ms,mt)
      T0_x  = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
      T0_y  = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
      T0_p  = eq_p(mp,6,ms,mt)
      T0_s  = eq_s(mp,6,ms,mt)
      T0_t  = eq_t(mp,6,ms,mt)
      T0_ss = eq_ss(mp,6,ms,mt)
      T0_tt = eq_tt(mp,6,ms,mt)
      T0_st = eq_st(mp,6,ms,mt)

      Vpar0    = eq_g(mp,7,ms,mt)
      Vpar0_x  = (   y_t(ms,mt) * eq_s(mp,7,ms,mt) - y_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
      Vpar0_y  = ( - x_t(ms,mt) * eq_s(mp,7,ms,mt) + x_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
      Vpar0_p  = eq_p(mp,7,ms,mt)
      Vpar0_s  = eq_s(mp,7,ms,mt)
      Vpar0_t  = eq_t(mp,7,ms,mt)
      Vpar0_ss = eq_ss(mp,7,ms,mt)
      Vpar0_st = eq_st(mp,7,ms,mt)
      Vpar0_tt = eq_tt(mp,7,ms,mt)

      P0    = abs(r0 * T0)
      P0_x  = r0_x * T0 + r0 * T0_x
      P0_y  = r0_y * T0 + r0 * T0_y
      P0_s  = r0_s * T0 + r0 * T0_s
      P0_t  = r0_t * T0 + r0 * T0_t
      P0_p  = r0_p * T0 + r0 * T0_p
      P0_ss = r0_ss * T0 + 2.d0 * r0_s * T0_s + r0 * T0_ss
      P0_tt = r0_tt * T0 + 2.d0 * r0_t * T0_t + r0 * T0_tt
      P0_st = r0_st * T0 + r0_s * T0_t + r0_t * T0_s + r0 * T0_st

      ps0_xx = (ps0_ss * y_t(ms,mt)**2 - 2.d0*ps0_st * y_s(ms,mt)*y_t(ms,mt) + ps0_tt * y_s(ms,mt)**2 &
             + ps0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                             &
             + ps0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2              &
             - xjac_x * (ps0_s* y_t(ms,mt) - ps0_t * y_s(ms,mt))  / xjac**2

      ps0_yy = (ps0_ss * x_t(ms,mt)**2 - 2.d0*ps0_st * x_s(ms,mt)*x_t(ms,mt) + ps0_tt * x_s(ms,mt)**2 &
             + ps0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                             &
             + ps0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2              &
             - xjac_y * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )  / xjac**2

      ps0_xy = (- ps0_ss * y_t(ms,mt)*x_t(ms,mt) - ps0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
              + ps0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - ps0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - ps0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2              &
              - xjac_x * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )   / xjac**2

      u0_xx = (u0_ss * y_t(ms,mt)**2 - 2.d0*u0_st * y_s(ms,mt)*y_t(ms,mt) + u0_tt * y_s(ms,mt)**2     &
            + u0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                               &
            + u0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )      / xjac**2              &
            - xjac_x * (u0_s * y_t(ms,mt) - u0_t * y_s(ms,mt)) / xjac**2

      u0_yy = (u0_ss * x_t(ms,mt)**2 - 2.d0*u0_st * x_s(ms,mt)*x_t(ms,mt) + u0_tt * x_s(ms,mt)**2     &
            + u0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                               &
            + u0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2              &
            - xjac_y * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) ) / xjac**2

      u0_xy = (- u0_ss * y_t(ms,mt)*x_t(ms,mt) - u0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + u0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                            &
              - u0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                            &
              - u0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2              &
              - xjac_x * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) )   / xjac**2
 
      w0_xy = (- w0_ss * y_t(ms,mt)*x_t(ms,mt) - w0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + w0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                            &
              - w0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                            &
              - w0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2               &
              - xjac_x * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )   / xjac**2
              
      w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2     &
             + w0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
             + w0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )     / xjac**2              &
             - xjac_x * (w0_s* y_t(ms,mt) - w0_t * y_s(ms,mt))  / xjac**2

      w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2     &
             + w0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
             + w0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2              &
             - xjac_y * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )  / xjac**2

      w0_2xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2 ) / xjac**2
      w0_2yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2 ) / xjac**2

      r0_xx = (r0_ss * y_t(ms,mt)**2 - 2.d0*r0_st * y_s(ms,mt)*y_t(ms,mt) + r0_tt * y_s(ms,mt)**2     &
            + r0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                               &
            + r0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2                &
            - xjac_x * (r0_s* y_t(ms,mt) - r0_t * y_s(ms,mt))  / xjac**2

      r0_yy = (r0_ss * x_t(ms,mt)**2 - 2.d0*r0_st * x_s(ms,mt)*x_t(ms,mt) + r0_tt * x_s(ms,mt)**2     &
            + r0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                               &
            + r0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2             &
            - xjac_y * (- r0_s * x_t(ms,mt) + r0_t * x_s(ms,mt) )  / xjac**2

      r0_xy = (- r0_ss * y_t(ms,mt)*x_t(ms,mt) - r0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + r0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                            &
              - r0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                            &
              - r0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2               &
              - xjac_x * (- r0_s * x_t(ms,mt) + r0_t * x_s(ms,mt) )   / xjac**2

      T0_xx = (T0_ss * y_t(ms,mt)**2 - 2.d0*T0_st * y_s(ms,mt)*y_t(ms,mt) + T0_tt * y_s(ms,mt)**2     &
            + T0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                               &
            + T0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2                &
            - xjac_x * (T0_s * y_t(ms,mt) - T0_t * y_s(ms,mt))  / xjac**2

      T0_yy = (T0_ss * x_t(ms,mt)**2 - 2.d0*T0_st * x_s(ms,mt)*x_t(ms,mt) + T0_tt * x_s(ms,mt)**2     &
            + T0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                               &
            + T0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2             &
            - xjac_y * (- T0_s * x_t(ms,mt) + T0_t * x_s(ms,mt) )  / xjac**2

      T0_xy = (- T0_ss * y_t(ms,mt)*x_t(ms,mt) - T0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + T0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                            &
              - T0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                            &
              - T0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2              &
              - xjac_x * (- T0_s * x_t(ms,mt) + T0_t * x_s(ms,mt) )   / xjac**2

      vpar0_xx = (vpar0_ss * y_t(ms,mt)**2 - 2.d0*vpar0_st * y_s(ms,mt)*y_t(ms,mt) + vpar0_tt * y_s(ms,mt)**2  &
               + vpar0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                  &
               + vpar0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )   / xjac**2                    &
               - xjac_x * (vpar0_s * y_t(ms,mt) - vpar0_t * y_s(ms,mt)) / xjac**2

      vpar0_yy = (vpar0_ss * x_t(ms,mt)**2 - 2.d0*vpar0_st * x_s(ms,mt)*x_t(ms,mt) + vpar0_tt * x_s(ms,mt)**2  &
               + vpar0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                                  &
               + vpar0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2                 &
               - xjac_y * (- vpar0_s * x_t(ms,mt) + vpar0_t * x_s(ms,mt) ) / xjac**2

      vpar0_xy = (- vpar0_ss * y_t(ms,mt)*x_t(ms,mt) - vpar0_tt * x_s(ms,mt)*y_s(ms,mt)               &
                 + vpar0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                      &
                 - vpar0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                      &
                 - vpar0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2        &
                 - xjac_x * (- vpar0_s * x_t(ms,mt) + vpar0_t * x_s(ms,mt) )   / xjac**2

      P0_xx = r0_xx * T0 + 2.d0 * r0_x * T0_x + r0 * T0_xx
      P0_yy = r0_yy * T0 + 2.d0 * r0_y * T0_y + r0 * T0_yy
      P0_xy = r0_xy * T0 + r0_x * T0_y + r0_y * T0_x + r0 * T0_xy
      
      delta_u    = delta_g(mp,2,ms,mt)
      delta_u_s  = delta_s(mp,2,ms,mt)
      delta_u_t  = delta_t(mp,2,ms,mt)
      delta_u_ss = delta_ss(mp,2,ms,mt)
      delta_u_st = delta_st(mp,2,ms,mt)
      delta_u_tt = delta_tt(mp,2,ms,mt)

      delta_u_x = (   y_t(ms,mt) * delta_s(mp,2,ms,mt) - y_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac
      delta_u_y = ( - x_t(ms,mt) * delta_s(mp,2,ms,mt) + x_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac

      delta_ps_x = (   y_t(ms,mt) * delta_s(mp,1,ms,mt) - y_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
      delta_ps_y = ( - x_t(ms,mt) * delta_s(mp,1,ms,mt) + x_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
      
      delta_u_xx = (delta_u_ss * y_t(ms,mt)**2 - 2.d0*delta_u_st * y_s(ms,mt)*y_t(ms,mt)           &
                 + delta_u_tt * y_s(ms,mt)**2                                                      &
                 + delta_u_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                  &
                 + delta_u_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2     &
                 - xjac_x * (delta_u_s* y_t(ms,mt) - delta_u_t * y_s(ms,mt))         / xjac**2

      delta_u_yy = (delta_u_ss * x_t(ms,mt)**2 - 2.d0*delta_u_st * x_s(ms,mt)*x_t(ms,mt)           &
                 + delta_u_tt * x_s(ms,mt)**2                                                      &
                 + delta_u_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                  &
                 + delta_u_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )  / xjac**2     &
                 - xjac_y * (- delta_u_s * x_t(ms,mt) + delta_u_t * x_s(ms,mt) )     / xjac**2
     
      ! --- Temperature dependent resistivity
      if ( eta_T_dependent ) then
        eta_T     = eta   * (corr_neg_temp(T0)/T_0)**(-1.5d0)
        deta_dT   = - eta   * (1.5d0)  * corr_neg_temp(T0)**(-2.5d0) * T_0**(1.5d0)
        d2eta_d2T =   eta   * (3.75d0) * corr_neg_temp(T0)**(-3.5d0) * T_0**(1.5d0)
      else
        eta_T     = eta
        deta_dT   = 0.d0
        d2eta_d2T = 0.d0
      end if
     
      ! --- Temperature dependent viscosity
      if ( visco_T_dependent ) then
        visco_T   = visco * (corr_neg_temp(T0)/T_0)**(-1.5d0)
        dvisco_dT = - visco * (1.5d0)  * corr_neg_temp(T0)**(-2.5d0) * T_0**(1.5d0)
      else
        visco_T   = visco
        dvisco_dT = 0.d0
      end if
     
      ! --- Temperature dependent parallel heat diffusivity
      if ( ZKpar_T_dependent ) then
        ZKpar_T   = ZK_par * (corr_neg_temp(T0)/T_0)**(+2.5d0)        ! temperature dependent parallel conductivity
        dZKpar_dT = ZK_par * (2.5d0)  * corr_neg_temp(T0)**(+1.5d0) * T_0**(-2.5d0)
        if (ZKpar_T .gt. ZK_par_max) then
          ZKpar_T   = Zk_par_max
          dZKpar_dT = 0.d0
        endif
      else
        ZKpar_T   = ZK_par                                            ! parallel conductivity
        dZKpar_dT = 0.d0
      endif

      eta_num_T   = eta_num                                           ! hyperresistivity
      visco_num_T = visco_num                                         ! hyperviscosity



          ! --- Vorticity Equation Terms calculation (with 2 at the end means strong form 
          ! --- without the 2 at the end is the weak form expression) 
          ! --- The terms calculated in these 2 different ways should give the same result, verify!

          rot_NL_Ve      = - r0_hat * BigR**2. * w0 * (v_s(mp,ms,mt) * u0_t - v_t(mp,ms,mt) * u0_s)
          
          rot_NL_Ve_2    = - v(mp,ms,mt) * ( BigR**4. * r0 * ( u0_x * w0_y - u0_y * w0_x )         &
                                           - BigR**3. * 2.* r0 * w0 * u0_y )                  * xjac

          gd_rho_NL_Ve   = - 0.5d0 * vv2 * (v_x(mp,ms,mt) * r0_y_hat                               &
                           - v_y(mp,ms,mt) * r0_x_hat)                                        * xjac
                          
          gd_rho_NL_Ve_2 = - v(mp,ms,mt) * ( BigR**2. * r0_x_hat * ( u0_yy * u0_y + u0_xy * u0_x ) &
                                               - BigR * r0_y_hat * (  u0_y * u0_y + u0_x  * u0_x   &
                                               + BigR * ( u0_xy * u0_y + u0_xx * u0_x)  )          &
                                    + BigR**2. * w0 * ( u0_x * r0_y_hat - u0_y * r0_x_hat ) ) * xjac

          Poiss_Psi_j    = + v(mp,ms,mt) * ( ps0_s * zj0_t - ps0_t * zj0_s )
          
          Poiss_Psi_j_2  = + v(mp,ms,mt) * ( ps0_x * zj0_y - ps0_y * zj0_x )                  * xjac
          
          der_j_tor      = - v(mp,ms,mt) * eps_cyl * F0 / BigR * zj0_p                        * xjac
          
          rot_grad_P     = + BigR**2. * ( v_s(mp,ms,mt) * p0_t - v_t(mp,ms,mt) * p0_s )
          
          rot_grad_P_2   = - BigR * 2. * v(mp,ms,mt) * p0_y                                   * xjac
          
          visco_lap_W    = - BigR * visco_T  * (v_x(mp,ms,mt) * w0_x + v_y(mp,ms,mt) * w0_y)  * xjac
          
          visco_lap_W_2  = + v(mp,ms,mt) * BigR * ( visco_T * ( w0_xx + w0_x/BigR + w0_yy ) ) * xjac
          
          visco_num_W    = - visco_num_T * (v_xx(mp,ms,mt) + v_x(mp,ms,mt)/Bigr + v_yy(mp,ms,mt))  &
                           *(w0_xx + w0_x/BigR + w0_yy)                                       * xjac
                           
          visco_num_W_2  = - visco_num_T * (v_xx(mp,ms,mt) + v_x(mp,ms,mt)/Bigr + v_yy(mp,ms,mt))  &
                           *(w0_2xx + w0_x/BigR + w0_2yy)                                     * xjac
          
          rot_Vstar      = - v(mp,ms,mt) * tauIC * BigR**4. * (p0_s * w0_t - p0_t * w0_s)          &
                           - tauIC * BigR**3. * p0_y * (v_x(mp,ms,mt) * u0_x                       &
                                                      + v_y(mp,ms,mt) * u0_y)              * xjac  &
                           - v(mp,ms,mt) * tauIC * BigR**4. * (u0_xy * (p0_xx - p0_yy)             &
                                                             - p0_xy * (u0_xx - u0_yy) )   * xjac
!                          -3*v*tauIC* BigR**3 * (p0_x*u0_xy - p0_y*u0_xx) * xjac

          rot_Vstar_2    = - v(mp,ms,mt) * ( tauIC*BigR**4. *( p0_x * w0_y - p0_y * w0_x )         &
                           + tauIC*BigR**4. *( u0_xy * (p0_xx - p0_yy) - p0_xy * (u0_xx - u0_yy) ) &
                           - tauIC*BigR**3. *( p0_xy*u0_x + p0_y*u0_xx + p0_yy*u0_y + p0_y*u0_yy ) &
                           - tauIC*BigR**2. *( p0_y * u0_x ) )                                * xjac
!                      -3*v*tauIC* BigR**2 * (p0_x*u0_xy - p0_y*u0_xx) * xjac * tstep &

          rot_delta_u    = -1. * ( BigR * r0_hat * (   v_x(mp,ms,mt) * delta_u_x                   &
                                               + v_y(mp,ms,mt) * delta_u_y ) * xjac ) / tstep_n(1)
                                             
          rot_delta_u_2  = v(mp,ms,mt) * (   BigR**3.     * r0   * delta_u_xx                      &
                                           + BigR**3.     * r0_x * delta_u_x                       &
                                           + 3.* BigR**2. * r0   * delta_u_x                       &
                                           + BigR**3.     * r0   * delta_u_yy                      &
                                           + BigR**3.     * r0_y * delta_u_y   )                   &
                                           / tstep_n(1)                                       * xjac

          Source_term    = + BigR**3. * particle_source(ms,mt) * (v_x(mp,ms,mt) * u0_x             &
                                                                + v_y(mp,ms,mt) * u0_y)       * xjac
                                                                
          psi_norm = get_psi_n( ps0, y_g(ms,mt)) 
        
          Source_term_2  = -1. * v(mp,ms,mt) * ( particle_source(ms,mt) * BigR**3. * w0            &
                           + BigR * ( u0_x * ps0_x + u0_y * ps0_y )                                &
                           / ( sqrt( ps0_x * ps0_x + ps0_y * ps0_y ) )                             &
                           * ( particle_source(ms,mt) * 2. * BigR * ps0_x                          &
                           / ( sqrt( ps0_x * ps0_x + ps0_y * ps0_y ) )                             &
                           + BigR**2. * ( 0.5 * particlesource *                                   &
                           (tanh( (psi_norm - particlesource_psin) / particlesource_sig ) - 1.)))) &
                           * xjac
          
          
          psi_abs  = sqrt(ps0_x * ps0_x + ps0_y * ps0_y)
          Btheta   = (psi_abs / BigR)
          Btot     = sqrt(F0**2. + ps0_x**2. + ps0_y**2.) / BigR
          
          ! --- Storage of the integration of the Terms in one element and average values in the element
          RHS_local(1,1,1)  = RHS_local(1,1,1)  + BigR
          RHS_local(2,1,1)  = RHS_local(2,1,1)  + BigZ
          RHS_local(3,1,1)  = RHS_local(3,1,1)  + ps0
          RHS_local(4,1,1)  = RHS_local(4,1,1)  + u0
          RHS_local(5,1,1)  = RHS_local(5,1,1)  + zj0
          RHS_local(6,1,1)  = RHS_local(6,1,1)  + w0
          RHS_local(7,1,1)  = RHS_local(7,1,1)  + r0
          RHS_local(8,1,1)  = RHS_local(8,1,1)  + T0
          RHS_local(9,1,1)  = RHS_local(9,1,1)  + Vpar0
          RHS_local(10,1,1) = RHS_local(10,1,1) - (u0_x * ps0_x + u0_y * ps0_y) / psi_abs
          
          RHS_local(11,1,1) = RHS_local(11,1,1) + -1./Btheta * ( (u0_x                             &
                                                        + tauIC/r0*(T0_x*r0 + r0_x*T0)) * ps0_x    &
                                                        +        (u0_y                             &
                                                        + tauIC/r0*(T0_y*r0 + r0_y*T0)) * ps0_y )  &
                                                        + Vpar0 * Btheta
                                                        
          RHS_local(12,1,1) = RHS_local(12,1,1) + sqrt(GAMMA*T0) / Btot
          RHS_local(13,1,1) = RHS_local(13,1,1) + Btot
          RHS_local(14,1,1) = RHS_local(14,1,1) + psi_norm 
          
          RHS_local(15,1,1) = RHS_local(15,1,1) + -1./Btheta * ( (u0_x * ps0_x) + (u0_y * ps0_y) )
          
          RHS_local(16,1,1) = RHS_local(16,1,1) + -1./Btheta * ( (tauIC/r0*(T0_x*r0 + r0_x*T0) * ps0_x)  &
                                                               + (tauIC/r0*(T0_y*r0 + r0_y*T0) * ps0_y))
                                                           
          RHS_local(17,1,1) = RHS_local(17,1,1) + Vpar0 * Btheta
          
          RHS_local(18,1,1) = RHS_local(18,1,1) + rot_delta_u    * wst
          RHS_local(19,1,1) = RHS_local(19,1,1) + rot_delta_u_2  * wst
          RHS_local(20,1,1) = RHS_local(20,1,1) + Poiss_Psi_j    * wst
          RHS_local(21,1,1) = RHS_local(21,1,1) + Poiss_Psi_j_2  * wst
          RHS_local(22,1,1) = RHS_local(22,1,1) + rot_grad_P     * wst
          RHS_local(23,1,1) = RHS_local(23,1,1) + rot_grad_P_2   * wst
          RHS_local(24,1,1) = RHS_local(24,1,1) + der_j_tor      * wst
          RHS_local(25,1,1) = RHS_local(25,1,1) + rot_NL_Ve      * wst
          RHS_local(26,1,1) = RHS_local(26,1,1) + rot_NL_Ve_2    * wst
          RHS_local(27,1,1) = RHS_local(27,1,1) + gd_rho_NL_Ve   * wst
          RHS_local(28,1,1) = RHS_local(28,1,1) + gd_rho_NL_Ve_2 * wst
          RHS_local(29,1,1) = RHS_local(29,1,1) + visco_lap_W    * wst
          RHS_local(30,1,1) = RHS_local(30,1,1) + visco_lap_W_2  * wst
          RHS_local(31,1,1) = RHS_local(31,1,1) + visco_num_W    * wst
          RHS_local(32,1,1) = RHS_local(32,1,1) + visco_num_W_2  * wst
          RHS_local(33,1,1) = RHS_local(33,1,1) + Source_term    * wst
          RHS_local(34,1,1) = RHS_local(34,1,1) + Source_term_2  * wst
          RHS_local(35,1,1) = RHS_local(35,1,1) + v_s(mp,ms,mt)  * wst
          RHS_local(36,1,1) = RHS_local(36,1,1) + v_t(mp,ms,mt)  * wst
          RHS_local(37,1,1) = RHS_local(37,1,1) + xjac           * wst
          RHS_local(38,1,1) = RHS_local(38,1,1) + rot_Vstar      * wst
          RHS_local(39,1,1) = RHS_local(39,1,1) + rot_Vstar_2    * wst
          RHS_local(40,1,1) = RHS_local(40,1,1) + r0 * vv2       * wst * BigR * xjac
          RHS_local(41,1,1) = RHS_local(41,1,1) + ( Poiss_Psi_j + rot_grad_P + der_j_tor ) * wst
          RHS_local(42,1,1) = RHS_local(42,1,1) + (   Poiss_Psi_j + rot_grad_P + der_j_tor         &
                                                    + visco_lap_W + visco_num_W_2        ) * wst
                                                    
          RHS_local(43,1,1) = RHS_local(43,1,1) + particle_source(ms,mt)
          
    enddo
  enddo
 enddo
 
 ! --- Grid for VTK
 xyz(1:3,ife) = (/ ( RHS_local(1,1,1) / float(n_gauss*n_gauss*n_plane_local) )                     &
                 , ( RHS_local(2,1,1) / float(n_gauss*n_gauss*n_plane_local) ), 0.d0 /)
 
 ! --- Storage of the integration of the terms in one element for VTK
 scalars(ife,1)  = RHS_local(3,1,1)  / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,2)  = RHS_local(4,1,1)  / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,3)  = RHS_local(5,1,1)  / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,4)  = RHS_local(6,1,1)  / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,5)  = RHS_local(7,1,1)  / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,6)  = RHS_local(8,1,1)  / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,7)  = RHS_local(9,1,1)  / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,8)  = RHS_local(10,1,1) / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,9)  = RHS_local(11,1,1) / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,10) = RHS_local(12,1,1) / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,11) = RHS_local(13,1,1) / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,12) = RHS_local(14,1,1) / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,13) = RHS_local(15,1,1) / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,14) = RHS_local(16,1,1) / float(n_gauss*n_gauss*n_plane_local)
 scalars(ife,15) = RHS_local(17,1,1) / float(n_gauss*n_gauss*n_plane_local)
 
 scalars(ife,16) = RHS_local(18,1,1)
 scalars(ife,17) = RHS_local(20,1,1)
 scalars(ife,18) = RHS_local(22,1,1)
 scalars(ife,19) = RHS_local(24,1,1)
 scalars(ife,20) = RHS_local(25,1,1)
 scalars(ife,21) = RHS_local(27,1,1)
 scalars(ife,22) = RHS_local(29,1,1)
 scalars(ife,23) = RHS_local(32,1,1)
 scalars(ife,24) = RHS_local(33,1,1)
 scalars(ife,25) = RHS_local(38,1,1)
 scalars(ife,26) = RHS_local(40,1,1)
 scalars(ife,27) = RHS_local(41,1,1)
 scalars(ife,28) = RHS_local(42,1,1)
 
 scalars(ife,29) = RHS_local(43,1,1) / float(n_gauss*n_gauss*n_plane_local)
 

 ! --- Sum of the integrated values (of each vorticity term)
 ! --- for all the elements respecting the following conditions (if not commented)
  
!  if ( ( ife .GT. (n_tht + 2) )  .AND.            &
!       ( ife .LT. (n_flux * n_tht - 4 * n_tht) ) )  then

            Sum_rot_delta_u  = Sum_rot_delta_u  + RHS_local(18,1,1)
            Sum_Poiss_Psi_j  = Sum_Poiss_Psi_j  + RHS_local(20,1,1)
            Sum_rot_grad_P   = Sum_rot_grad_P   + RHS_local(22,1,1)
            Sum_rot_grad_P_2 = Sum_rot_grad_P_2 + RHS_local(23,1,1)
            Sum_der_j_tor    = Sum_der_j_tor    + RHS_local(24,1,1)
            Sum_rot_NL_Ve    = Sum_rot_NL_Ve    + RHS_local(25,1,1)
            Sum_gd_rho_NL_Ve = Sum_gd_rho_NL_Ve + RHS_local(27,1,1)
            Sum_visco_lap_W  = Sum_visco_lap_W  + RHS_local(29,1,1)
            Sum_visco_num_W  = Sum_visco_num_W  + RHS_local(32,1,1)
            Sum_Source_term  = Sum_Source_term  + RHS_local(33,1,1)
            Sum_rot_Vstar    = Sum_rot_Vstar    + RHS_local(38,1,1)
            
            Sum_xjac   = Sum_xjac  + RHS_local(37,1,1)
            Sum_vv2    = Sum_vv2   + RHS_local(40,1,1)

            Sum_sum_1  = Sum_sum_1 + RHS_local(41,1,1)
            Sum_sum_2  = Sum_sum_2 + RHS_local(42,1,1)
            
            count_Ve = count_Ve + 1
            
            RHS_local(44,1,1) = 1.

!endif
 
 write(999992, '(1I7,29ES16.7)')  ife               &
                                , RHS_local(1,1,1)  &
                                , RHS_local(2,1,1)  &
                                , RHS_local(3,1,1)  &
                                , RHS_local(18,1,1) &
                                , RHS_local(19,1,1) &
                                , RHS_local(20,1,1) &
                                , RHS_local(21,1,1) &
                                , RHS_local(22,1,1) &
                                , RHS_local(23,1,1) &
                                , RHS_local(24,1,1) &
                                , RHS_local(25,1,1) &
                                , RHS_local(26,1,1) &
                                , RHS_local(27,1,1) &
                                , RHS_local(28,1,1) &
                                , RHS_local(29,1,1) &
                                , RHS_local(30,1,1) &
                                , RHS_local(31,1,1) &
                                , RHS_local(32,1,1) &
                                , RHS_local(33,1,1) &
                                , RHS_local(34,1,1) &
                                , RHS_local(35,1,1) &
                                , RHS_local(36,1,1) &
                                , RHS_local(37,1,1) &
                                , RHS_local(38,1,1) &
                                , RHS_local(39,1,1) &
                                , RHS_local(40,1,1) &
                                , RHS_local(41,1,1) &
                                , RHS_local(42,1,1) &
                                , RHS_local(44,1,1)

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   END LOOP ELEMENTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

close (unit = 999992)

print*, ' '
print*, 'count_Ve  = ', count_Ve
print*, 'count_Ve2 = ', count_Ve2
print*, ' '

! -- Possible extra normalization (not used if commented)
!Sum_rot_NL_Ve    = Sum_rot_NL_Ve    / float(count_Ve) / Sum_xjac
!Sum_gd_rho_NL_Ve = Sum_gd_rho_NL_Ve / float(count_Ve) / Sum_xjac
!Sum_Poiss_Psi_j  = Sum_Poiss_Psi_j  / float(count_Ve) / Sum_xjac
!Sum_der_j_tor    = Sum_der_j_tor    / float(count_Ve) / Sum_xjac
!Sum_rot_grad_P   = Sum_rot_grad_P   / float(count_Ve) / Sum_xjac
!Sum_rot_grad_P_2 = Sum_rot_grad_P_2 / float(count_Ve) / Sum_xjac
!Sum_visco_lap_W  = Sum_visco_lap_W  / float(count_Ve) / Sum_xjac
!Sum_visco_num_W  = Sum_visco_num_W  / float(count_Ve) / Sum_xjac
!Sum_rot_Vstar    = Sum_rot_Vstar    / float(count_Ve) / Sum_xjac
!Sum_Source_term  = Sum_Source_term  / float(count_Ve) / Sum_xjac
!Sum_rot_delta_u  = Sum_rot_delta_u  / float(count_Ve) / Sum_xjac

!Sum_vv2          = Sum_vv2          / float(count_Ve) / Sum_xjac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! --- WRITE the vorticity terms averaged in the considered group of elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! --- Check exist and Open file
      write(filename,'(a)')  'Average_Terms_Voriticity_Gauss.dat'
        
      ! --- Check existance
      inquire(file=filename, exist=exist_file)

      If  ( exist_file ) Then
      
        write(*,*) 'File Average_Terms_Voriticity_Gauss.dat exist'
        
        open(unit=file_handle, file=filename, status='old'             &
           , position='append', action='write', iostat=errorFile)
            
        if (errorFile /= 0) then
          write(*,*) 'ERROR : Opening file Average_Terms_Voriticity_Gauss"' &
                    , trim(filename), '" failed.'
          stop
        end if

                                       
        write(file_handle,'(15ES16.7)') t_start           &
                                      , Sum_rot_NL_Ve     &
                                      , Sum_gd_rho_NL_Ve  &
                                      , Sum_Poiss_Psi_j   &
                                      , Sum_der_j_tor     &
                                      , Sum_rot_grad_P    &
                                      , Sum_rot_grad_P_2  &
                                      , Sum_visco_lap_W   &
                                      , Sum_visco_num_W   &
                                      , Sum_rot_Vstar     &
                                      , Sum_Source_term   &
                                      , Sum_rot_delta_u   &
                                      , Sum_vv2           &
                                      , Sum_sum_1         &
                                      , Sum_sum_2
        
      ElseIf( .NOT. exist_file) Then
      
        write(*,*) 'File does not exist then create'

        open(unit=file_handle, file=filename, status='new'    &
            ,action='write', iostat=errorFile)
            
        if (errorFile /= 0) then
          write(*,*) 'ERROR : Creating file Average_Terms_Voriticity_Gauss"' &
                    , trim(filename), '" failed.'
          stop
        end if

        write(file_handle,'(15A16)') ' n2 t_start     '   &
                                   , ' n3 Ve_NL       '   &
                                   , ' n4 rho_NL      '   &
                                   , ' n5 PoissPsiJ   '   &
                                   , ' n6 der_j_tor   '   &
                                   , ' n7 grdP_weak   '   &
                                   , ' n8 grdP_strg   '   &
                                   , ' n9 visco_W     '   &
                                   , ' n10 visco_Num  '   &
                                   , ' n11 Vstar      '   &
                                   , ' n12 Src        '   &
                                   , ' n13 Sum_delta  '   &
                                   , ' n14 Sum_vv2    '   &
                                   , ' n15 SumEqui    '   &
                                   , ' n16 SumEquiVis ' 

                                       
        write(file_handle,'(15ES16.7)') t_start           &
                                      , Sum_rot_NL_Ve     &
                                      , Sum_gd_rho_NL_Ve  &
                                      , Sum_Poiss_Psi_j   &
                                      , Sum_der_j_tor     &
                                      , Sum_rot_grad_P    &
                                      , Sum_rot_grad_P_2  &
                                      , Sum_visco_lap_W   &
                                      , Sum_visco_num_W   &
                                      , Sum_rot_Vstar     &
                                      , Sum_Source_term   &
                                      , Sum_rot_delta_u   &
                                      , Sum_vv2           &
                                      , Sum_sum_1         &
                                      , Sum_sum_2
                                      
      Endif

      ! --- Close file
      close (unit=file_handle)
      
      write(*,*) 'END WRITE'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   END WRITE vorticity terms averaged in the considered group of elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! --- Write the binary VTK file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Find elements neighbours for VTK element definition
allocate(element_neighbours(3,element_list%n_elements))

element_neighbours = 0

n_element_less_private = n_flux * n_tht + n_open * (n_tht) !+ 2*(n_leg - 1))

print*, ' '
print*, ' n_element_less_private = ', n_element_less_private
print*, ' '

do i=1,(element_list%n_elements-1)

  if ( i .LE. n_element_less_private ) then 
  
    if ( ( (mod(i,n_tht) .EQ. 0) .OR. (mod((i-1),n_tht) .EQ. 0) .OR. (mod((i+1),n_tht) .EQ. 0) )   &
        .AND. (i .LE. (n_flux * n_tht)) ) then
  
      do j=(i-n_tht+1),(i+n_tht+1)
        if ( neighbours(node_list,element_list%element(i),element_list%element(j),iside_i,iside_j) ) then
          element_neighbours(iside_i,i) = j
        endif
      enddo
  
    else
  
      do j=(i+1),(i + n_tht + 2*n_leg + 4) !element_list%n_elements
        if ( neighbours(node_list,element_list%element(i),element_list%element(j),iside_i,iside_j) ) then
          element_neighbours(iside_i,i) = j
        endif
      enddo
    
    endif
    
  else
  
    do j=( i - n_open*(n_tht + 2*(n_leg)) ),element_list%n_elements
      if ( neighbours(node_list,element_list%element(i),element_list%element(j),iside_i,iside_j) ) then
        element_neighbours(iside_i,i) = j
      endif
    enddo
    
  endif
enddo

! --- Some prints for check
!print*, ' '
!print*, ' element_neighbours(:,1)             = ', element_neighbours(:,1)
!print*, ' element_neighbours(:,2)             = ', element_neighbours(:,2)
!print*, ' element_neighbours(:,200)           = ', element_neighbours(:,200)
!print*, ' element_neighbours(:,199)           = ', element_neighbours(:,199)
!print*, ' element_neighbours(:,201)           = ', element_neighbours(:,201)
!print*, ' element_neighbours(:,26217)         = ', element_neighbours(:,26217)
!print*, ' element_neighbours(:,26218)         = ', element_neighbours(:,26218)

!print*, ' '
!print*, ' element_list%element(1)%vertex(:)   = ', element_list%element(1)%vertex(:)
!print*, ' element_list%element(2)%vertex(:)   = ', element_list%element(2)%vertex(:)
!print*, ' element_list%element(199)%vertex(:) = ', element_list%element(199)%vertex(:)
!print*, ' element_list%element(200)%vertex(:) = ', element_list%element(200)%vertex(:)
!print*, ' element_list%element(201)%vertex(:) = ', element_list%element(201)%vertex(:)
!print*, ' '
!print*, ' 1 node_list%node(1)%x(1,1,1)          = ', node_list%node(element_list%element(1)%vertex(1))%x(1,1,1)
!print*, ' 1 node_list%node(2)%x(1,1,1)          = ', node_list%node(element_list%element(1)%vertex(2))%x(1,1,1)
!print*, ' 1 node_list%node(3)%x(1,1,1)          = ', node_list%node(element_list%element(1)%vertex(3))%x(1,1,1)
!print*, ' 1 node_list%node(4)%x(1,1,1)          = ', node_list%node(element_list%element(1)%vertex(4))%x(1,1,1)
!print*, ' '
!print*, ' 2 node_list%node(1)%x(1,1,1)          = ', node_list%node(element_list%element(2)%vertex(1))%x(1,1,1)
!print*, ' 2 node_list%node(2)%x(1,1,1)          = ', node_list%node(element_list%element(2)%vertex(2))%x(1,1,1)
!print*, ' 2 node_list%node(3)%x(1,1,1)          = ', node_list%node(element_list%element(2)%vertex(3))%x(1,1,1)
!print*, ' 2 node_list%node(4)%x(1,1,1)          = ', node_list%node(element_list%element(2)%vertex(4))%x(1,1,1)
!print*, ' '
!print*, ' 1 node_list%node(1)%x(1,1,2)          = ', node_list%node(element_list%element(1)%vertex(1))%x(1,1,2)
!print*, ' 1 node_list%node(2)%x(1,1,2)          = ', node_list%node(element_list%element(1)%vertex(2))%x(1,1,2)
!print*, ' 1 node_list%node(3)%x(1,1,2)          = ', node_list%node(element_list%element(1)%vertex(3))%x(1,1,2)
!print*, ' 1 node_list%node(4)%x(1,1,2)          = ', node_list%node(element_list%element(1)%vertex(4))%x(1,1,2)
!print*, ' '
!print*, ' 2 node_list%node(1)%x(1,1,2)          = ', node_list%node(element_list%element(2)%vertex(1))%x(1,1,2)
!print*, ' 2 node_list%node(2)%x(1,1,2)          = ', node_list%node(element_list%element(2)%vertex(2))%x(1,1,2)
!print*, ' 2 node_list%node(3)%x(1,1,2)          = ', node_list%node(element_list%element(2)%vertex(3))%x(1,1,2)
!print*, ' 2 node_list%node(4)%x(1,1,2)          = ', node_list%node(element_list%element(2)%vertex(4))%x(1,1,2)
!print*, ' '


do i=1,element_list%n_elements

  if ( (element_neighbours(1,i) .GT. 0) .AND. (element_neighbours(2,i) .GT. 0) .AND.               &
       (element_neighbours(3,i) .GT. 0) ) then
  
    ien(1,i) = i - 1
    ien(2,i) = element_neighbours(1,i) - 1
    ien(3,i) = element_neighbours(2,i) - 1
    ien(4,i) = element_neighbours(3,i) - 1
    
    count_Ve2 = count_Ve2 + 1
    
  endif
enddo

print*, ' '
print*, '2 count_Ve  = ', count_Ve
print*, '2 count_Ve2 = ', count_Ve2
print*, ' '

count_Ve2 = 0

do i=1,element_list%n_elements
  
  if ( (element_neighbours(1,i) .GT. 0) .AND. (element_neighbours(2,i) .GT. 0) .AND.               &
       (element_neighbours(3,i) .GT. 0) ) then
         
    count_Ve2 = count_Ve2 + 1

    ien2(1,count_Ve2) = ien(1,i)
    ien2(2,count_Ve2) = ien(2,i)
    ien2(3,count_Ve2) = ien(3,i)
    ien2(4,count_Ve2) = ien(4,i)
        
  endif
enddo

print*, ' '
print*, '3 count_Ve  =  ' , count_Ve
print*, '3 count_Ve2 = ' , count_Ve2
print*, ' '


etype = 9  ! for VTK_QUAD

lf = char(10) ! line feed character


#ifdef IBM_MACHINE
open(unit=ivtk,file='jorek_tmp.vtk',form='unformatted',access='stream',status='replace')
#else
open(unit=ivtk,file='jorek_tmp.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN',status='replace')
#endif


buffer = '# vtk DataFile Version 3.0'//lf    ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                    ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                        ; write(ivtk) trim(buffer)
buffer = 'DATASET UNSTRUCTURED_GRID'//lf     ; write(ivtk) trim(buffer)

! POINTS SECTION
write(str1(1:12),'(i12)') (element_list%n_elements)
buffer = 'POINTS '//str1//'  float'//lf      ; write(ivtk) trim(buffer)
write(ivtk) ((real(xyz(i,j),4),i=1,3),j=1,(element_list%n_elements))

! CELLS SECTION
write(str1(1:12),'(i12)') (count_Ve2)    ! number of elements (cells)
write(str2(1:12),'(i12)') 5*(count_Ve2)  ! size of the following element list
buffer = lf//'CELLS '//str1//' '//str2//lf   ; write(ivtk) trim(buffer)
write(ivtk) (int(4,4),(int(ien2(j,i),4),j=1,4),i=1,(count_Ve2))

! CELL_TYPES SECTION
write(str1(1:12),'(i12)') (count_Ve2)    ! number of elements (cells)
buffer = lf//'CELL_TYPES'//str1//lf          ; write(ivtk) trim(buffer)
write(ivtk) (int(etype,4),i=1,(count_Ve2))

! POINT_DATA SECTION
write(str1(1:12),'(i12)') (element_list%n_elements)
buffer = lf//'POINT_DATA '//str1             ; write(ivtk) trim(buffer)

do i_var = 1, n_scalars
  buffer = lf//'SCALARS '//scalar_names(i_var)//' float'//lf ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf
  write(ivtk) trim(buffer)
  write(ivtk) (real(scalars(i,i_var),4),i=1,(element_list%n_elements))
enddo

close(ivtk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! --- END write the binary VTK file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*) 'Done.'
 
end program jorek2vtk_GaussVortTerms
