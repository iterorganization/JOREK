!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2vtk_3d

use constants
use data_structure
use phys_module
use mod_import_restart
use mod_interp
implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
real*8,allocatable    :: HZ(:,:)
integer,allocatable   :: ien (:,:)
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, index, index_node, n_points
integer               :: n_toroidal
character             :: buffer*80, lf*1, str1*10, str2*10
character*8, allocatable :: scalar_names(:), vector_names(:)
real*4                :: float
real*8                :: s, t, phi, angle, cur_pert
real*8                :: P,P_s,P_t,P_st,P_ss,P_tt,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8                :: Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt, ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt, W,W_s,W_t,W_st,W_ss,W_tt
real*8                :: U,U_s,U_t,U_st,U_ss,U_tt, RHO,RH_s,RH_t,RH_st,RH_ss,RH_tt, TT,TT_s,TT_t,TT_st,TT_ss,TT_tt
real*8                :: u0_x, u0_y, xjac, v_perp, Psi_J, R_p, error, zj_x, zj_y, ps_x, ps_y
logical               :: periodic, density_only
integer               :: ierr, my_id
logical               :: without_n0_mode, RphiZ_coords

namelist /vtk_params/ nsub, without_n0_mode, periodic, RphiZ_coords

write(*,*) 'jorek2vtk_3d'

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub            = 5        		! Number of subdivisions of the cubic finite elements into linear pieces
without_n0_mode = .false.  		! If true, do not include the n=0 mode (i_tor=1)
periodic        = .true.		! Are we doing the whole tor?
density_only    = .false.		! Write density only (for smaller vtk file)
n_toroidal      = 200 !n_plane 		! Number of toroidal snapshots
RphiZ_coords    = .false.               ! use xyz transformation from JOREK wiki 

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

! --- Number of scalars and vectors written to the VTK file
if(density_only) then
  n_scalars = 1
  n_vectors = 0
else
  n_scalars = 7
  n_vectors = 1
endif

do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
enddo

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
nnos = n_toroidal * nsub*nsub*node_list%n_nodes

allocate(xyz(3,nnos), scalars(nnos,1:n_scalars), scalar_names(n_scalars))
if(density_only) then
  scalar_names = (/'density '/)
else
  allocate(vector_names(n_vectors), vectors(nnos,3,1:n_vectors))
  vector_names = (/ 'B_field '/)
  scalar_names = (/ 'flux    ','U       ','j       ','omega   ','density ','T       ','v_par   '/)
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

do m=1,n_toroidal
  if (periodic) then
    phi = 2.d0 * PI * float(m-1)/float(n_toroidal)
  else
    phi = 2.d0 * PI * float(m-1)/float(n_toroidal-1) / float(n_period)
  endif
  HZ(1,m)   = 1.d0
  do i=1,(n_tor-1)/2
    HZ(2*i,m)     = cos(mode(2*i)  *phi)
    HZ(2*i+1,m)   = sin(mode(2*i+1)*phi)
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

        call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
        xjac  = R_s * Z_t - R_t * Z_s
        if ( xjac == 0.d0 ) xjac = 1.d-8 ! (workaround to avoid floating invalid)

        inode = inode+1

        if (RphiZ_coords) then
          xyz(1:3,inode) = (/ R * cos(angle), -R*sin(angle), Z /)   !from the JOREK wiki
        else
          xyz(1:3,inode) = (/ R * cos(angle), Z, R*sin(angle) /)
        endif

        do i_tor = 1,n_tor

          if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode
         
          if(density_only) then
            call interp(node_list,element_list,i,var_rho,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,1) = scalars(inode,1) + P * HZ(i_tor,m)
	  else
            call interp(node_list,element_list,i,var_psi,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,1) = scalars(inode,1) + P * HZ(i_tor,m)

            ps_x  = (	Z_t * P_s - Z_s * P_t ) / xjac
            ps_y  = ( - R_t * P_s + R_s * P_t ) / xjac

            vectors(inode,1:3,1) = vectors(inode,1:3,1) + (/+ ps_y * HZ(i_tor,m) / R * cos(angle),	  &
            						    - ps_x * HZ(i_tor,m) / R,			  &
            						    + ps_y * HZ(i_tor,m) / R * sin(angle)  /)
            if (i_tor .eq. 1) then
              vectors(inode,1:3,1) = vectors(inode,1:3,1) + (/ - F0/R * sin(angle), -0.d0, F0/R *cos(angle)/)
            endif

            call interp(node_list,element_list,i,var_u,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,2) = scalars(inode,2) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_zj,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,3) = scalars(inode,3) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_w,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,4) = scalars(inode,4) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_rho,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,5) = scalars(inode,5) + P * HZ(i_tor,m)

            call interp(node_list,element_list,i,var_T,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,6) = scalars(inode,6) + P * HZ(i_tor,m)
            
            if (with_Vpar) then
              call interp(node_list,element_list,i,var_Vpar,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            else
              P = 0.d0
            end if
            scalars(inode,var_Vpar) = scalars(inode,var_Vpar) + P * HZ(i_tor,m)
	  endif

	enddo

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
