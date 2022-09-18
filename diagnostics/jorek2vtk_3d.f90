!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2vtk_3d

use mod_parameters, only: n_order
use mod_chi
use constants
use data_structure
use phys_module
use mod_import_restart
use mod_interp
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
integer               :: n_scalar_forces, n_scalars_base, n_vectors_base 
integer               :: n_forces_scalars, n_forces_vectors, n_gvec_scalars, n_gvec_vectors
integer               :: s_scalars_forces, s_vectors_forces, s_scalars_gvec, s_vectors_gvec
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
real*8,allocatable    :: HZ(:,:), HZ_p(:,:), HZ_pp(:,:), HZ_coord(:,:), HZ_coord_p(:,:)
integer,allocatable   :: ien (:,:)
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, index, index_node, n_points, k_tor
integer               :: n_toroidal
character             :: buffer*80, lf*1, str1*10, str2*10
character*8, allocatable :: scalar_names(:), vector_names(:)
real*4                :: float
real*8                :: s, t, phi, angle, cur_pert
real*8                :: P,P_s,P_t,P_st,P_ss,P_tt
real*8                :: R,R_s,R_t,R_phi,R_st,R_ss,R_tt,R_sp,R_tp,R_pp
real*8                :: Z,Z_s,Z_t,Z_p,Z_st,Z_ss,Z_tt,Z_sp,Z_tp,Z_pp
real*8                :: Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt, ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt, ZJ_phi, W,W_s,W_t,W_st,W_ss,W_tt
real*8                :: ps_p, ps_phi, ps_sx, ps_sy, ps_tx, ps_ty, ps_xx, ps_yy, ps_xy, ps_xp, ps_yp,ps_pp, ps_yphi, ps_xphi, ps_phiphi, ps_x_itor, ps_y_itor
real*8                :: Lap_ps
real*8                :: delta_phi
real*8                :: U,U_s,U_t,U_st,U_ss,U_tt, RHO,RH_s,RH_t,RH_st,RH_ss,RH_tt, TT,TT_s,TT_t,TT_st,TT_ss,TT_tt
real*8                :: u_x, u_y, u_p
real*8                :: rho_x, rho_y, rho_phi
real*8                :: TT_x, TT_y, TT_phi
real*8                :: u0_x, u0_y, xjac, xjac_s, xjac_t, xjac_x, xjac_y, v_perp, Psi_J, R_p, error, zj_x, zj_y, ps_x, ps_y
real*8                :: BR, BZ, BP
real*8                :: BR_gvec, BRg,BRg_s,BRg_t,BRg_st,BRg_ss,BRg_tt, BR_R, BR_Z, BR_P
real*8                :: BZ_gvec, BZg,BZg_s,BZg_t,BZg_st,BZg_ss,BZg_tt, BZ_R, BZ_Z, BZ_P
real*8                :: BP_gvec, Bpg,Bpg_s,Bpg_t,Bpg_st,Bpg_ss,Bpg_tt, BP_R, BP_Z, BP_P
real*8                :: JR_gvec, JZ_gvec, JP_gvec, PR_gvec, PZ_gvec, PP_gvec, FR_gvec, FZ_gvec, FP_gvec
real*8                :: PR, PZ, PP, JR, JZ, JP, FR, FZ, FP
real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

logical               :: periodic, density_only
integer               :: ierr, my_id
logical               :: without_n0_mode, RphiZ_coords, include_forces, include_gvec_fields

namelist /vtk_params/ nsub, without_n0_mode, periodic, RphiZ_coords, density_only, n_toroidal, include_forces, include_gvec_fields

write(*,*) 'jorek2vtk_3d'

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub            = 3        		! Number of subdivisions of the cubic finite elements into linear pieces
without_n0_mode = .false.  		! If true, do not include the n=0 mode (i_tor=1)
periodic        = .false.		! Are we doing the whole tor?
density_only    = .false.		! Write density only (for smaller vtk file)
n_toroidal      = 6 !n_plane 		! Number of toroidal snapshots
RphiZ_coords    = .false.               ! use xyz transformation from JOREK wiki 
include_forces  = .true.               ! calculate the force balance of the JOREK fields
include_gvec_fields    = .true.               ! include the imported fields from GVEC (only makes sense for stellarator simulations) 

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
write(*,*) 'without_n0_mode =', without_n0_mode
write(*,*) 'periodic        =', periodic
write(*,*)
if (include_gvec_fields) then
#if ((JOREK_MODEL == 183) || (JOREK_MODEL == 083))
  continue
#else
  write(*,*) 'GVEC fields can only be written for stellarator models!'
  stop
#endif
endif

! --- Number of scalars and vectors written to the VTK file
if(density_only) then
  n_scalars_base = 1
  n_vectors_base = 0
  n_scalars = n_scalars_base
  n_vectors = n_vectors_base
else
  n_scalars_base = 7
  n_vectors_base = 2
  n_scalars = n_scalars_base
  n_vectors = n_vectors_base
endif

if (include_forces) then
  s_scalars_forces = n_scalars
  n_forces_scalars = 1
  n_scalars        = n_scalars + n_forces_scalars
  s_vectors_forces = n_vectors
  n_forces_vectors = 3
  n_vectors = n_vectors+n_forces_vectors
endif

if (include_gvec_fields) then
  s_scalars_gvec = n_scalars
  s_vectors_gvec = n_vectors
  n_gvec_scalars = 1
  n_gvec_vectors = 4
  n_scalars = n_scalars+n_gvec_scalars
  n_vectors = n_vectors+n_gvec_vectors
endif

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
enddo
                                                       
do k_tor=1, n_coord_tor
  mode_coord(k_tor) = + int(k_tor / 2) * n_coord_period
enddo

call init_chi_basis


call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
nnos = n_toroidal * nsub*nsub*node_list%n_nodes

allocate(xyz(3,nnos), scalars(nnos,1:n_scalars), scalar_names(n_scalars))
if(density_only) then
  scalar_names(1:n_scalars_base) = (/'density '/)
else
  allocate(vector_names(n_vectors), vectors(nnos,3,1:n_vectors))
  vector_names(1:n_vectors_base) = (/ 'B_field' , 'v_field'/)
  scalar_names(1:n_scalars_base) = (/ 'flux    ','U       ','j       ','omega   ','density ','T       ','F      '/)
endif

if (include_forces) then
  scalar_names(s_scalars_forces+1:s_scalars_forces+n_forces_scalars) = (/ 'delta_phi  '/)
  vector_names(s_vectors_forces+1:s_vectors_forces+n_forces_vectors) = (/ 'J_field', 'gradP_field', 'F_field'/)
endif
if (include_gvec_fields) then
  vector_names(s_vectors_gvec+1:s_vectors_gvec+n_gvec_vectors) = (/ 'B_gvec', 'J_gvec', 'gradP_gvec', 'F_gvec'/)
  scalar_names(s_scalars_gvec+1:s_scalars_gvec+n_gvec_scalars) = (/ 'j_gvec' /)
endif

if (periodic) then
  nel   = (n_toroidal)   * (nsub-1)*(nsub-1)*element_list%n_elements
else
  nel   = (n_toroidal-1) * (nsub-1)*(nsub-1)*element_list%n_elements
endif

nnoel = 8
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0
n_points = nsub*nsub*element_list%n_elements        ! number of points in one poloidal plane

allocate(HZ(n_tor,n_toroidal))
allocate(HZ_p(n_tor,n_toroidal))
allocate(HZ_pp(n_tor,n_toroidal))
allocate(HZ_coord(n_coord_tor,n_toroidal))
allocate(HZ_coord_p(n_coord_tor,n_toroidal))

do m=1,n_toroidal
  if (periodic) then
    phi = 2.d0 * PI * float(m-1)/float(n_toroidal)
  else
    phi = 2.d0 * PI * float(m-1)/float(n_toroidal-1) / float(n_period)
  endif
  HZ(1,m)   = 1.d0
  do i=1,(n_tor-1)/2
    HZ(2*i,m)     =                        cos(mode(2*i)  *phi)
    HZ_p(2*i,m)   = -float(mode(2*i))    * sin(mode(2*i)  *phi)
    HZ_pp(2*i,m)  = -float(mode(2*i))**2 * cos(mode(2*i)  *phi)
    HZ(2*i+1,m)   =                        sin(mode(2*i+1)*phi)
    HZ_p(2*i+1,m) =  float(mode(2*i))    * cos(mode(2*i+1)*phi)
    HZ_pp(2*i+1,m)= -float(mode(2*i))**2 * sin(mode(2*i+1)*phi)
  enddo

  HZ_coord(1,m) = 1.0
  HZ_coord_p(1,k) = 0.d0
  do i=1,(n_coord_tor-1)/2
    HZ_coord(2*i,m)      =                           cos(mode_coord(2*i)  *phi)
    HZ_coord_p(2*i,m)    = - float(mode_coord(2*i))      * sin(mode_coord(2*i)  *phi)
    HZ_coord(2*i+1,m)    =                         - sin(mode_coord(2*i+1)*phi)
    HZ_coord_p(2*i+1,m)  = - float(mode_coord(2*i+1))    * cos(mode_coord(2*i+1)*phi)
  enddo
enddo

do m=1, n_toroidal
  ! --- Print progress information as jorek2vtk_3d may run very long...
  if ( mod(m,n_toroidal/40+1) == 0 ) write(*,'(" Plane ",i4.4," of ",i4.4)') m, n_toroidal

  if (periodic) then
    angle = 2.d0 * PI * float(m-1)/float(n_toroidal)
  else
    angle = 2.d0 * PI * float(m-1)/float(n_toroidal-1) / float(n_period)
  endif

  do i=1,element_list%n_elements

    do j=1,nsub
      s = float(j-1)/float(nsub-1)
      do k=1,nsub
        t = float(k-1)/float(nsub-1)

        ! The following 50 lines could be replaced with interp_PRZ(_1) (after adding without_n0_mode there, or manually subtracting)

        call interp_RZP(node_list,element_list,i,s,t,angle,R,R_s,R_t,R_phi,R_st,R_ss,R_tt,R_sp,R_tp,R_pp, &
                       Z,Z_s,Z_t,Z_p,Z_st,Z_ss,Z_tt,Z_sp,Z_tp,Z_pp)
        chi = get_chi(R,Z,angle)
        inode = inode+1
         
        xjac  = R_s * Z_t - R_t * Z_s
        xjac_s = R_ss*Z_t + R_s*Z_st - R_st*Z_s - R_t*Z_ss
        xjac_t = R_st*Z_t + R_s*Z_tt - R_tt*Z_s - R_t*Z_st
        xjac_x = (xjac_s*Z_t - xjac_t*Z_s) / xjac
        xjac_y = (R_s*xjac_t - R_t*xjac_s) / xjac
        if ( xjac == 0.d0 ) xjac = 1.d-8 ! (workaround to avoid floating invalid)

        if (RphiZ_coords) then
          xyz(1:3,inode) = (/ R * cos(angle), -R*sin(angle), Z /)   !from the JOREK wiki
        else
          xyz(1:3,inode) = (/ R * cos(angle), Z, R*sin(angle) /)
        endif

        ps_x  = 0.d0; ps_y  = 0.d0; ps_p  = 0.d0; ps_xx = 0.d0; ps_yy = 0.d0; ps_xy = 0.d0; ps_yp = 0.d0; ps_xp = 0.d0; ps_pp = 0.d0
        zj    = 0.d0; zj_y  = 0.d0; zj_x  = 0.d0; zj_phi  = 0.d0
        do i_tor = 1,n_tor

          if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode
          if (n_coord_period .ne. 1 .and. without_n0_mode .and. mod(mode(i_tor),n_coord_period) .eq. 0) cycle
         
          if(density_only) then
            call interp(node_list,element_list,i,var_rho,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,1) = scalars(inode,1) + P * HZ(i_tor,m)
          else
            call interp(node_list,element_list,i,var_psi,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,1) = scalars(inode,1) + P * HZ(i_tor,m)

            ps_x  = ps_x + (  Z_t * P_s - Z_s * P_t  ) / xjac * HZ(i_tor,m)
            ps_y  = ps_y + ( -R_t * P_s + R_s * P_t ) / xjac * HZ(i_tor,m)
            if ((jorek_model .eq. 083) .or. (jorek_model .eq. 183)) then
                ! Get psi_p for B_JOREK field representation
                ps_p = ps_p + P*HZ_p(i_tor,m) - (   Z_t * P_s - Z_s * P_t ) / xjac * HZ(i_tor,m)*R_phi &
                                              - ( - R_t * P_s + R_s * P_t ) / xjac * HZ(i_tor,m)*Z_p

                ! Get higher derivatives for J_JOREK representation
                ps_sx = ((Z_st * P_s + Z_t * P_ss - Z_ss * P_t - Z_s * P_st) * xjac - xjac_s * (  Z_t * P_s - Z_s * P_t  )) / xjac ** 2
                ps_tx = ((Z_tt * P_s + Z_t * P_st - Z_st * P_t - Z_s * P_tt) * xjac - xjac_t * (  Z_t * P_s - Z_s * P_t  )) / xjac ** 2 
                ps_xx = ps_xx + ( Z_t * ps_sx - Z_s * ps_tx) / xjac * HZ(i_tor, m)
                
                ps_sy = (( -R_st * P_s - R_t * P_ss + R_ss * P_t + R_s * P_st) * xjac - xjac_s * ( -R_t * P_s + R_s * P_t  )) / xjac ** 2
                ps_ty = (( -R_tt * P_s - R_t * P_st + R_st * P_t + R_s * P_tt) * xjac - xjac_t * ( -R_t * P_s + R_s * P_t  )) / xjac ** 2 
                ps_yy = ps_yy + ( -R_t * ps_sy + R_s * ps_ty) / xjac * HZ(i_tor, m)
                ps_xy = ps_xy + ( -R_t * ps_sx + R_s * ps_tx ) / xjac * HZ(i_tor,m)
                
                ps_xp = ps_xp + (  Z_t * P_s - Z_s * P_t  ) / xjac * HZ_p(i_tor,m)
                ps_yp = ps_yp + ( -R_t * P_s + R_s * P_t  ) / xjac * HZ_p(i_tor,m)
                ps_pp = ps_pp +  P * HZ_pp(i_tor,m)
            endif

            call interp(node_list,element_list,i,var_u,i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
            scalars(inode,2) = scalars(inode,2) + U * HZ(i_tor,m)

            u_x  = u_x   + (   Z_t * U_s - Z_s * U_t )     / xjac * HZ(i_tor,m)
            u_y  = u_y   + ( - R_t * U_s + R_s * U_t )     / xjac * HZ(i_tor,m)
            if ((jorek_model .eq. 083) .or. (jorek_model .eq. 183)) then
                u_p = u_p + U*HZ_p(i_tor,m) - (   Z_t * U_s - Z_s * U_t ) / xjac * HZ(i_tor,m)*R_phi &
                                            - ( - R_t * U_s + R_s * U_t ) / xjac * HZ(i_tor,m)*Z_p
            endif

            call interp(node_list,element_list,i,var_zj,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,3) = scalars(inode,3) + P * HZ(i_tor,m)
            zj   = zj     +  P * HZ(i_tor,m)   
            zj_x = zj_x   + (   Z_t * P_s - Z_s * P_t )     / xjac * HZ(i_tor,m)  
            zj_y = zj_y   + ( - R_t * P_s + R_s * P_t )     / xjac * HZ(i_tor,m)
            zj_phi = zj_phi   + P * HZ_p(i_tor,m) - R_phi * (   Z_t * P_s - Z_s * P_t )     / xjac * HZ(i_tor,m) &
                                              - Z_p   * ( - R_t * P_s + R_s * P_t )     / xjac * HZ(i_tor,m)


            call interp(node_list,element_list,i,var_w,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,4) = scalars(inode,4) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_rho,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,5) = scalars(inode,5) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_T,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,6) = scalars(inode,6) + P * HZ(i_tor,m)
            
            call interp(node_list,element_list,i,var_F,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,7) = scalars(inode,var_F) + P * HZ(i_tor,m)
          endif
        enddo
        
        if (     (jorek_model .eq. 083).or. (jorek_model .eq. 183) )then
            BR = chi(1,0,0)      + (ps_y*chi(0,0,1) - ps_p*chi(0,1,0))/(F0*R)
            BZ = chi(0,1,0)      - (ps_x*chi(0,0,1) - ps_p*chi(1,0,0))/(F0*R)
            BP = chi(0,0,1)/R    + (ps_x*chi(0,1,0) - ps_y*chi(1,0,0))/F0       
            vectors(inode,1:3, 1) = (/ BR * cos(angle) - BP * sin(angle), &
                                       BZ, &
                                       BR * sin(angle) + BP * cos(angle) /)

            vectors(inode,1:3, 2) = (/ -R * u_y * cos(angle) - u_p * sin(angle), &
                                        R * u_x, &
                                       -R * u_y * sin(angle) + u_p * cos(angle) /)
        else
            vectors(inode,1:3, 1) = (/- F0/R * sin(angle) + ps_y / R * cos(angle),       &
                                                          - ps_x / R,                     &
                                        F0/R * cos(angle) + ps_y / R * sin(angle)/)

            vectors(inode,1:3, 2) = (/                      -R * u_y * cos(angle) , &
                                                             R * u_x, &
                                                            -R * u_y * sin(angle)  /)
        endif
        
        if (include_gvec_fields) then 
          ! Interpolate force balance from GVEC
          BR_gvec = 0; BZ_gvec = 0; BP_gvec = 0
          JR_gvec = 0; JZ_gvec = 0; JP_gvec = 0
          do i_tor=1, n_coord_tor
           call interp_gvec(node_list,element_list,i,1,1,i_tor,s,t,BRg,BRg_s,BRg_t,BRg_st,BRg_ss,BRg_tt)
           call interp_gvec(node_list,element_list,i,1,2,i_tor,s,t,BZg,BZg_s,BZg_t,BZg_st,BZg_ss,BZg_tt)
           call interp_gvec(node_list,element_list,i,1,3,i_tor,s,t,Bpg,Bpg_s,Bpg_t,Bpg_st,Bpg_ss,Bpg_tt)
           BR_gvec  = BR_gvec + BRg * HZ_coord(i_tor, m)          
           BZ_gvec  = BZ_gvec + BZg * HZ_coord(i_tor, m)          
           BP_gvec  = BP_gvec + BPg * HZ_coord(i_tor, m)          
          
           BR_R  = (   Z_t * BRg_s - Z_s * BRg_t )     / xjac * HZ_coord(i_tor,m)
           BR_Z  = ( - R_t * BRg_s + R_s * BRg_t )     / xjac * HZ_coord(i_tor,m)
           BR_p  = BRg * HZ_coord_p(i_tor, m) - BR_R * R_phi - Z_p * BR_Z
           BZ_R  = (   Z_t * BZg_s - Z_s * BZg_t )     / xjac * HZ_coord(i_tor,m)
           BZ_Z  = ( - R_t * BZg_s + R_s * BZg_t )     / xjac * HZ_coord(i_tor,m)
           BZ_p  = BZg * HZ_coord_p(i_tor, m) - BZ_R * R_phi - Z_p * BZ_Z
           BP_R  = (   Z_t * BPg_s - Z_s * BPg_t )     / xjac * HZ_coord(i_tor,m)
           BP_Z  = ( - R_t * BPg_s + R_s * BPg_t )     / xjac * HZ_coord(i_tor,m)
           BP_p  = BPg * HZ_coord_p(i_tor, m) - BP_R * R_phi - Z_p * BP_Z
           
           JR_gvec = JR_gvec + (1 / R * BZ_p - BP_Z)
           JZ_gvec = JZ_gvec + (1 / R * (BPg + R * BP_R - BR_p))
           JP_gvec = JP_gvec + (BR_Z - BZ_R)
          enddo
          scalars(inode, s_scalars_gvec+1) = F0 / (chi(1,0,0)**2 + chi(0,1,0)**2 + chi(0,0,1)**2/R**2) * (chi(1,0,0) * JR_gvec + chi(0,1,0) * JZ_gvec + chi(0,0,1)/R * JP_gvec)
          vectors(inode,:, s_vectors_gvec+1) = (/ BR_gvec * cos(angle) - BP_gvec * sin(angle), &
                                   BZ_gvec, &
                                   BR_gvec * sin(angle) + BP_gvec * cos(angle) /)
          
          
          vectors(inode,:, s_vectors_gvec+2) = (/ JR_gvec * cos(angle) - JP_gvec * sin(angle), &
                                   JZ_gvec, &
                                   JR_gvec * sin(angle) + JP_gvec * cos(angle) /)
          
          call interp_gvec(node_list,element_list,i,3,1,i_tor,s,t,BRg,BRg_s,BRg_t,BRg_st,BRg_ss,BRg_tt)
          PR_gvec = BRg_s * Z_t / xjac
          PZ_gvec = BRg_s * (-1) * R_t / xjac
          PP_gvec = BRg_s * (R_t * Z_p - Z_t * R_phi)
          vectors(inode,:, s_vectors_gvec+3) = (/ PR_gvec * cos(angle) - PP_gvec * sin(angle), &
                                   PZ_gvec, &
                                   PR_gvec * sin(angle) + PP_gvec * cos(angle) /)        

          FR_gvec = (JP_gvec * BZ_gvec - JZ_gvec * BP_gvec) / MU_ZERO  - PR_gvec 
          FZ_gvec = (-JR_gvec * BP_gvec + JP_gvec * BR_gvec) / MU_ZERO - PZ_gvec
          FP_gvec = (JZ_gvec * BR_gvec - JR_gvec * BZ_gvec) / MU_ZERO  - PP_gvec 
          vectors(inode,:, s_vectors_gvec+4) = (/ FR_gvec * cos(angle) - FP_gvec * sin(angle), &
                                   FZ_gvec, &
                                   FR_gvec * sin(angle) + FP_gvec * cos(angle) /)
        endif

        if (include_forces) then
#if ((JOREK_MODEL == 183) || (JOREK_MODEL == 083))
          ! Get Jorek magnetic field
          BR = chi(1,0,0)      + (ps_y*chi(0,0,1) - ps_p*chi(0,1,0))/(F0*R)
          BZ = chi(0,1,0)      - (ps_x*chi(0,0,1) - ps_p*chi(1,0,0))/(F0*R)
          BP = chi(0,0,1)/R    + (ps_x*chi(0,1,0) - ps_y*chi(1,0,0))/F0
          
          ! Get Jorek current
          ps_phi    = ps_p
          ps_xphi   = ps_xp - R_phi*ps_xx - Z_p*ps_xy
          ps_yphi   = ps_yp - R_phi*ps_xy - Z_p*ps_yy
          ps_phiphi = ps_pp - R_phi * ps_xp - Z_p * ps_yp 
          
          Lap_ps    = ps_xx + ps_x/R + ps_yy + ps_phiphi/R**2

          ! j = (grad(chi).grad)grad(psi)    - Lap(psi)     -     (grad(psi).grad)grad(chi)
          JR = (chi(1,0,0)*ps_xx + chi(0,0,1)/R**2*ps_xphi + chi(0,1,0)*ps_xy                            - chi(1,0,0) * Lap_ps   &
               - (chi(2,0,0)*ps_x + chi(1,0,1)/R**2*ps_phi + chi(1,1,0)*ps_y)) / F0         
          JZ = (chi(1,0,0)*ps_xy + chi(0,0,1)/R**2*ps_yphi + chi(0,1,0)*ps_yy                            - chi(0,1,0) * Lap_ps   &
               - (chi(1,1,0)*ps_x + chi(0,1,1)/R**2*ps_phi + chi(0,2,0)*ps_y)) / F0
          JP = (chi(1,0,0)*(ps_xphi/R - ps_phi/R**2) + chi(0,0,1)*ps_phiphi/R**3 + chi(0,1,0)*ps_yphi/R  - chi(0,0,1) * Lap_ps/R &
               - (ps_x*(chi(1,0,1)/R - chi(0,0,1)/R**2) + chi(0,0,2)*ps_phi/R**3 + ps_y*chi(0,2,0))) / F0
          vectors(inode,:,s_vectors_forces+1) = (/  JR * cos(angle) - JP * sin(angle), &
                                                    JZ, &
                                                    JR * sin(angle) + JP * cos(angle) /)

          rho    = 0.0; rho_x   = 0.0;  rho_y  = 0.0;   rho_phi = 0.0
          TT     = 0.0;   TT_x  = 0.0;    TT_y = 0.0;   TT_phi  = 0.0
          do i_tor = 1,n_tor
            call interp(node_list,element_list,i,var_rho,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            rho      = rho     + P * HZ(i_tor,m) 
            rho_x    = rho_x   + (   Z_t * P_s - Z_s * P_t )     / xjac * HZ(i_tor,m)
            rho_y    = rho_y   + ( - R_t * P_s + R_s * P_t )     / xjac * HZ(i_tor,m)
            rho_phi  = rho_phi + P * HZ_p(i_tor,m) - R_phi * (   Z_t * P_s - Z_s * P_t )     / xjac * HZ(i_tor,m) &
                                                 - Z_p   * ( - R_t * P_s + R_s * P_t )     / xjac * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_T,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            TT      = TT     + P * HZ(i_tor,m) 
            TT_x    = TT_x   + (   Z_t * P_s - Z_s * P_t )     / xjac * HZ(i_tor,m)
            TT_y    = TT_y   + ( - R_t * P_s + R_s * P_t )     / xjac * HZ(i_tor,m)
            TT_phi  = TT_phi + P * HZ_p(i_tor,m) - R_phi * (   Z_t * P_s - Z_s * P_t )     / xjac * HZ(i_tor,m) &
                                               - Z_p   * ( - R_t * P_s + R_s * P_t )     / xjac * HZ(i_tor,m)
          enddo
          PR = rho_x * TT + rho * TT_x
          PZ = rho_y * TT + rho * TT_y
          PP = 1.0 / R * (rho_phi * TT + rho * TT_phi)
          vectors(inode,:,s_vectors_forces+2) = (/  PR * cos(angle) - PP * sin(angle), &
                                                    PZ, &
                                                    PR * sin(angle) + PP * cos(angle) /) 

          FR = ( JP*BZ - JZ*BP - PR) / MU_ZERO
          FZ = (-JR*BP + JP*BR - PZ) / MU_ZERO
          FP = ( JZ*BR - JR*BZ - PP) / MU_ZERO
          vectors(inode,:,s_vectors_forces+3) = (/ FR * cos(angle) - FP * sin(angle), &
                                                   FZ, &
                                                   FR * sin(angle) + FP * cos(angle) /)

          ! - v Bv_parderiv(zj0)
          delta_phi =           - (zj_x * chi(1,0,0) + zj_y * chi(0,1,0) + zj_phi * chi(0,0,1) / (R*R)) / F0
          ! - v*Bv_pbrack(zj0,Psi0)
          delta_phi = delta_phi - ((zj_y*ps_phi - zj_phi*ps_y)*chi(1,0,0) + (zj_phi*ps_x - zj_x*ps_phi)*chi(0,1,0) + (zj_x*ps_y - zj_y*ps_x)*chi(0,0,1)) / R / (F0*F0)
          ! + Bv_pbrack(v,rho0*T0))/Bv2
          !delta_phi = delta_phi + ((v_y*PP/R - v_p*PZ/R)*chi(1,0,0) + (v_p*PR/R - v_x*PP/R)*chi(0,1,0) + (v_x*PZ - v_y*PR)*chi(0,0,1) / R)/Bv2
          scalars(inode, s_scalars_forces+1) = delta_phi

#else       
#if fullmhd
          write(*,*) 'Inclusion of the JOREK forces is not implemented yet for full MHD!'
          stop
#endif

          ps_x   = 0.d0; ps_y   = 0.d0; ps_xp  = 0.d0;  ps_yp   = 0.d0; zj = 0.d0
          rho    = 0.0; rho_x   = 0.0;  rho_y  = 0.0;   rho_phi = 0.0
          TT     = 0.0;   TT_x  = 0.0;    TT_y = 0.0;   TT_phi  = 0.0
          do i_tor = 1,n_tor
            
            call interp(node_list,element_list,i,var_psi,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            
            ps_x   = ps_x + (  Z_t * P_s - Z_s * P_t  ) / xjac * HZ(i_tor,m)
            ps_y   = ps_y + ( - R_t * P_s + R_s * P_t ) / xjac * HZ(i_tor,m)
            ps_xp  = ps_xp + (  Z_t * P_s - Z_s * P_t  ) / xjac * HZ(i_tor,m) * HZ_p(i_tor, m)
            ps_yp  = ps_yp + ( - R_t * P_s + R_s * P_t ) / xjac * HZ(i_tor,m) * HZ_p(i_tor, m)
            
            call interp(node_list,element_list,i,var_zj,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            zj     = zj + P * HZ(i_tor, m)

            call interp(node_list,element_list,i,var_rho,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            rho    = rho     + P * HZ(i_tor,m) 
            rho_x  = rho_x   + (   Z_t * P_s - Z_s * P_t )     / xjac * HZ(i_tor,m)
            rho_y  = rho_y   + ( - R_t * P_s + R_s * P_t )     / xjac * HZ(i_tor,m)
            rho_phi  = rho_phi   + P * HZ_p(i_tor,m)

            call interp(node_list,element_list,i,var_T,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            TT    = TT     + P * HZ(i_tor,m) 
            TT_x  = TT_x   + (   Z_t * P_s - Z_s * P_t )     / xjac * HZ(i_tor,m)
            TT_y  = TT_y   + ( - R_t * P_s + R_s * P_t )     / xjac * HZ(i_tor,m)
            TT_phi  = TT_phi   + P * HZ_p(i_tor,m)
          enddo
          
          JR = 1.0 / R**2 * ps_xp
          JZ = 1.0 / R**2 * ps_yp
          JP = -1.0 / R * zj
          vectors(inode,:,s_vectors_forces+1) = (/  JR * cos(angle) - JP * sin(angle), &
                                                    JZ, &
                                                    JR * sin(angle) + JP * cos(angle) /) 

          PR = rho_x * TT + rho * TT_x
          PZ = rho_y * TT + rho * TT_y
          PP = 1.0 / R * (rho_phi * TT + rho * TT_phi)
          vectors(inode,:,s_vectors_forces+2) = (/  PR * cos(angle) - PP * sin(angle), &
                                                    PZ, &
                                                    PR * sin(angle) + PP * cos(angle) /) 

          FR = (-1.0 / R**2 * zj * ps_x - F0 / R**3 * ps_yp - PR) / MU_ZERO
          FZ = (F0 / R**3 * ps_xp - 1 / R**2 * zj * ps_y    - PZ) / MU_ZERO
          FP = (-1.0 / R**3 * (ps_x * ps_xp + ps_y*ps_yp)   - PP) / MU_ZERO
          vectors(inode,:,s_vectors_forces+3) = (/ FR * cos(angle) - FP * sin(angle), &
                                                   FZ, &
                                                   FR * sin(angle) + FP * cos(angle) /)
#endif
        endif

      enddo
    enddo

    if (m .lt. n_toroidal) then

      do j=1,nsub-1
        do k=1,nsub-1

          ielm        = ielm+1
          ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! 0 based indices for VTK
          ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
          ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
          ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k

          ien(5,ielm) = ien(1,ielm) + n_points
          ien(6,ielm) = ien(2,ielm) + n_points
          ien(7,ielm) = ien(3,ielm) + n_points
          ien(8,ielm) = ien(4,ielm) + n_points

        enddo
      enddo

    endif

    if ( (periodic) .and. (m .eq. n_toroidal)) then

      do j=1,nsub-1
        do k=1,nsub-1

          ielm        = ielm+1
          ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! 0 based indices for VTK
          ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
          ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
          ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k

          ien(5,ielm) = ien(1,ielm) - n_points * (n_toroidal-1)
          ien(6,ielm) = ien(2,ielm) - n_points * (n_toroidal-1)
          ien(7,ielm) = ien(3,ielm) - n_points * (n_toroidal-1)
          ien(8,ielm) = ien(4,ielm) - n_points * (n_toroidal-1)

        enddo
      enddo

    endif

  enddo
enddo


!--------------------------------------------------- write the binary VTK file
etype = 12  ! for vtk_quad

lf = char(10) ! line feed character

#ifdef IBM_MACHINE
open(unit=ivtk,file='jorek_tmp.vtk',form='unformatted',access='stream',status='replace')
#else
open(unit=ivtk,file='jorek_tmp.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN',status='replace')
#endif

buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)

! POINTS SECTION
write(str1(1:10),'(i10)') nnos
buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
write(ivtk) ((xyz(i,j),i=1,3),j=1,nnos)

! CELLS SECTION
write(str1(1:10),'(i10)') nel            ! number of elements (cells)
write(str2(1:10),'(i10)') nel*(1+nnoel)  ! size of the following element list (nel*(nnoel+1))
buffer = lf//lf//'CELLS '//str1//' '//str2//lf                                        ; write(ivtk) trim(buffer)
write(ivtk) (nnoel,(ien(i,j),i=1,nnoel),j=1,nel)

! CELL_TYPES SECTION
write(str1(1:10),'(i10)') nel   ! number of elements (cells)
buffer = lf//lf//'CELL_TYPES'//str1//lf                                               ; write(ivtk) trim(buffer)
write(ivtk) (etype,i=1,nel)

! POINT_DATA SECTION
write(str1(1:10),'(i10)') nnos
buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)

do i_var =1, n_scalars
  buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
  buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
  write(ivtk) (scalars(i,i_var),i=1,nnos)
enddo

if(.not. density_only) then
  do i_var =1, n_vectors
    buffer = lf//lf//'VECTORS '//vector_names(i_var)//' float'//lf                    ; write(ivtk) trim(buffer)
    write(ivtk) ((vectors(j,i,i_var),i=1,3),j=1,nnos)
  enddo
endif

close(ivtk)

write(*,*) 'done.'

end program jorek2vtk_3d
