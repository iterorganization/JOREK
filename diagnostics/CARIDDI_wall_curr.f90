!> Program to convert CARIDDI currents into vtk and povray
! Several files from the CARIDDI run are necessary:
! ndofel.dat
! gmat.dat
! iglobdof.dat
! S_mat.bin   
! x.dat
! ix.dat
! Then jorek_restart.h5 is converted into the vtk and povray
! Different components can be plotted separately (hardcoded)
! The correct mesh has to be chosen (nodes_per_elem = 6 or 8)
!===========================================================
program convert

  use hdf5
  use hdf5_io_module
  use constants
  
implicit none
integer, parameter :: n_components_max = 100

integer, parameter :: ind_max =    1000000
integer, parameter :: n_nodes_max = 100000
integer, parameter :: nodes_per_elem = 8
logical, parameter :: ascii = .false.
integer, parameter :: povray_colormap = 1
real*8,  parameter :: povray_wallcurr_min = 0.d0
real*8,  parameter :: povray_wallcurr_max = 1.d7
real*8             :: phimin = 2*pi/2., phimax = 2*pi
real*8  :: xyznode(n_nodes_max,3), xyzav(3), maxdist
integer :: elemnode(n_nodes_max/2, nodes_per_elem)
integer :: tmp_elemnode(nodes_per_elem)
integer :: elem_component(n_nodes_max/2)
integer :: n_nodes, n_elems, ierr, i, j, ifile = 143, min_elem, max_elem

integer :: nelems_comp, istart, iend
integer :: ncomponents, comp_start(n_components_max), comp_end(n_components_max)
character(len=80) :: buffer, dir_struct
character(len=3) :: ic
character(len=12) :: str1, str2
character(len=20) :: comp_name(n_components_max)
character(len=1), parameter :: lf = char(10)
character(len=20), parameter :: node_file='x.dat', elem_file='ix.dat'

logical :: include_struct(n_components_max) = .true.

! for current conversion
integer(HID_T)     :: file_id
integer :: n_wall_curr,l1,l2,ind_last, error
integer :: ij_glob(ind_max,2)=0, ind_glob(ind_max)=0, dof
real*8  :: xyz_e(3), phi
real*8, allocatable, dimension(:)     :: wall_curr, wall_curr_real, wall_curr_el
real*8, allocatable, dimension(:)     :: jx, jy, jz, jphi
real*8, allocatable, dimension(:,:)   :: S, gmat_sp
real*8, allocatable, dimension(:,:,:) :: gmat

integer, allocatable :: ndofel(:), ind_g(:,:), full_glob(:,:)


namelist /CARIDDI_plot/ phimin, phimax, ncomponents, comp_start, comp_end, comp_name, dir_struct

dir_struct = './'

ncomponents  = 1
comp_start   = 0
comp_end     = 0
comp_name(1) ='CARIDDI_all'

open(42, file='CARIDDI_plot.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from CARIDDI_plot.nml namelist.'
  read(42,CARIDDI_plot)
  close(42)
end if


min_elem = 1
max_elem = 100000

! --- Read nodes
open(42, file=trim(dir_struct)//trim(node_file), status='old', action='read')
i = 0
do
  i = i + 1
  read(42,*, iostat=ierr) xyznode(i,1), xyznode(i,2), xyznode(i,3)
  if ( ierr /= 0 ) exit
  write(98,'(3ES18.9)') xyznode(i,1), xyznode(i,3), xyznode(i,2) !###
end do
close(42)
n_nodes = i-1
write(*,*) 'Read ', n_nodes, ' nodes.'
elem_component=0
! --- Read elements
open(42, file=trim(dir_struct)//trim(elem_file), status='old', action='read')
i = 0
do
  i = i + 1
  read(42,*, iostat=ierr) elemnode(i,:), elem_component(i)
  if ( ierr /= 0 ) exit
  if ( nodes_per_elem == 8 ) then
    tmp_elemnode(:) = elemnode(i,:)
    elemnode(i,1) = tmp_elemnode(1)
    elemnode(i,2) = tmp_elemnode(5)
    elemnode(i,3) = tmp_elemnode(8)
    elemnode(i,4) = tmp_elemnode(4)
    elemnode(i,5) = tmp_elemnode(2)
    elemnode(i,6) = tmp_elemnode(6)
    elemnode(i,7) = tmp_elemnode(7)
    elemnode(i,8) = tmp_elemnode(3)
  elseif ( nodes_per_elem == 6 ) then
    tmp_elemnode(:) = elemnode(i,:)
    elemnode(i,1) = tmp_elemnode(1)
    elemnode(i,2) = tmp_elemnode(2)
    elemnode(i,3) = tmp_elemnode(3)
    elemnode(i,4) = tmp_elemnode(4)
    elemnode(i,5) = tmp_elemnode(5)
    elemnode(i,6) = tmp_elemnode(6)
  else
    write(*,*) 'STOP, not implemented'
    stop
  end if
  
  write(99,'(99i10)') elemnode(i,:), elem_component(i) !###
end do
close(42)
n_elems = i-1
write(*,*) 'Read ', n_elems, ' elements.'

! --- Check elements
write(*,*) '  Min node index in elements: ', minval(elemnode(1:n_elems,:))
write(*,*) '  Max node index in elements: ', maxval(elemnode(1:n_elems,:))
write(*,*) '  Indices first element: ', elemnode(1,:)
write(*,*) '  Indices last  element: ', elemnode(n_elems,:)
write(*,*)


!#============== LOAD number of DOFS per element
open(13, file=trim(dir_struct)//'ndofel.dat', action='read')
allocate(ndofel(n_elems))
do i=1,n_elems
  read(13,'(i12)') ndofel(i)
end do
close(13)

!=== gives current at baricenter of each element, i_element , i_edge, jx, jy, jz
open(13, file=trim(dir_struct)//'gmat.dat', action='read')
allocate(gmat_sp(ind_max,3),ind_g(ind_max,2))
ind_g=0; gmat_sp=0.d0
i=0
do
  i = i + 1
  read(13,*,iostat=ierr) ind_g(i,:), gmat_sp(i,:)
  if ( ierr /= 0 ) exit
  ind_last = i
end do
close(13)
! Transform into full matrix
l1=maxval(ind_g(:,1));l2=maxval(ind_g(:,2))
allocate(gmat(l1,l2,3))
gmat=0.d0
do i = 1,ind_last
  gmat(ind_g(i,1),ind_g(i,2),:) = gmat_sp(i,:)
end do
deallocate(gmat_sp)

! relationship between dof and element
open(13, file=trim(dir_struct)//'iglobdof.dat', action='read')
i=0
do
  i = i + 1
  read(13,*,iostat=ierr) ij_glob(i,:), ind_glob(i)
  if ( ierr /= 0 ) exit
  ind_last=i
end do
close(13)
l1=maxval(ij_glob(:,1));l2=maxval(ij_glob(:,2))

allocate(full_glob(l1,l2))
full_glob=0.d0
do i =1,ind_last
  full_glob(ij_glob(i,1),ij_glob(i,2)) = ind_glob(i)
end do

! ================= Convert wall currents =================================
call HDF5_open('jorek_restart.h5',file_id,error)
call HDF5_integer_reading(file_id,n_wall_curr,"n_wall_curr")

allocate(wall_curr(n_wall_curr),S(n_wall_curr,n_wall_curr),wall_curr_real(n_wall_curr))
call HDF5_array1D_reading(file_id,wall_curr,"wall_curr")
call HDF5_close(file_id)

! Load S matrix
open(13, file=trim(dir_struct)//'S_mat.bin', status='old', form='unformatted', access='stream', action='read')
read(13) S
close(13)
do i = 1, n_wall_curr
   wall_curr_real(i) = sum(S(i,:)*wall_curr)
end do
deallocate(wall_curr)

allocate(wall_curr_el(n_elems),jx(n_elems),jy(n_elems),jz(n_elems), jphi(n_elems))
wall_curr_el=0.d0

! ----------------- Calculate current density in element center of mass ----------------------------------
jx=0.d0;   jy=0.d0;    jz=0.d0
do i = 1, n_elems
  do j = 1, ndofel(i)
    dof = full_glob(i,j)
    jx(i) =jx(i) + gmat(i,j,1)*wall_curr_real(dof)/4e-7/PI
    jy(i) =jy(i) + gmat(i,j,2)*wall_curr_real(dof)/4e-7/PI
    jz(i) =jz(i) + gmat(i,j,3)*wall_curr_real(dof)/4e-7/PI
  end do
  do j = 1,nodes_per_elem
    xyz_e(:) = xyznode(elemnode(i,j),:)/nodes_per_elem
  end do
  phi = atan2(xyz_e(2),xyz_e(1))
  jphi(i) = -jx(i) * sin(phi) + jy(i) * cos(phi)
  wall_curr_el(i) = (jx(i)**2+jy(i)**2+jz(i)**2)**.5
end do



! --- export to VTK
!components to plot (0=all)
do i = 1, ncomponents
  if (( comp_start(i) .eq. -1) .and. (comp_end(i) .eq. -1)) then
    do j = 1, maxval(elem_component)
      write(ic,'(I0.3)') j
      call write_header('CARIDDI_comp'//trim(ic)//'.vtk',ifile, n_nodes, xyznode)
      call det_comp(elem_component, j, j,&
          n_elems, istart, iend, nelems_comp)
      call write_cell(ifile, n_elems, elemnode, istart, iend, nelems_comp)
      call write_header_scalar(ifile, nelems_comp)
      call write_scalar('j_w(Am^-3)', ifile, istart, iend, wall_curr_el)
      call write_scalar('I*e_phi', ifile, istart, iend, jphi)
      call write_vector('J', ifile, istart, iend, jx, jy, jz)
      close(ifile)
    end do
  else
    call write_header(trim(comp_name(i))//'.vtk',ifile, n_nodes, xyznode)
    call det_comp(elem_component, comp_start(i), comp_end(i),&
        n_elems, istart, iend, nelems_comp)
    call write_cell(ifile, n_elems, elemnode, istart, iend, nelems_comp)
    call write_header_scalar(ifile, nelems_comp)
    call write_scalar('j_w(Am^-3)', ifile, istart, iend, wall_curr_el)
    call write_scalar('I*e_phi', ifile, istart, iend, jphi)
    call write_vector('J', ifile, istart, iend, jx, jy, jz)
    close(ifile)
  end if
end do

! --- export to Povray
include_struct(:) = .false.
include_struct(1:6) = .true.
open(ifile, file='3dwall.pov', form='formatted', status='replace')
call write_header_pov(ifile)
call write_currdens_pov(ifile)
close(ifile)

contains
  subroutine write_header(filename, ifile, nnodes, xyznode)
    !> Write header for vtk file
    integer, intent(in)             :: nnodes, ifile
    real*8,  intent(in)             :: xyznode(:,:)
    character(len=*), intent(in)    :: filename
    character(len=80)               :: buffer

    integer           :: i, j
    integer           :: type
    character(len=12) :: str1, str2

    
    if ( ascii ) then
      open(ifile, file=filename, form='formatted', status='replace')
    else
      open(ifile, file=filename, form='unformatted', convert='BIG_ENDIAN', access='stream', status='replace')
    end if
    
    if ( ascii ) then
      write(ifile,'(a)')'# vtk DataFile Version 3.0'
      write(ifile,'(a)')'vtk output'
      write(ifile,'(a)')'ASCII'
      write(ifile,'(a)')'DATASET UNSTRUCTURED_GRID'
      write(str1,'(i12)') nnodes
      write(ifile,'(a)')'POINTS '//str1//'  float'
    else
      buffer = '# vtk DataFile Version 3.0'//lf    ; write(ifile) trim(buffer)
      buffer = 'vtk output'//lf                    ; write(ifile) trim(buffer)
      buffer = 'BINARY'//lf                        ; write(ifile) trim(buffer)
      buffer = 'DATASET UNSTRUCTURED_GRID'//lf     ; write(ifile) trim(buffer)
      write(str1,'(i12)') nnodes
      buffer = 'POINTS '//str1//'  float'//lf      ; write(ifile) trim(buffer)
    end if
    
    if ( ascii ) then
        write(ifile,'(99f12.5)') ((xyznode(j,i),i=1,3),j=1,nnodes)
    else
      write(ifile) ((real(xyznode(j,i),4),i=1,3),j=1,nnodes)
    end if
    
  end subroutine write_header

  subroutine write_cell(ifile, nelems, elemnode, istart, iend, nelems_final)
    !> write information between dof and node into file
    integer, intent(in)             :: nelems, ifile
    integer, intent(in)             :: istart, iend, nelems_final
    integer, intent(in)             :: elemnode(:,:)
    character(len=80)               :: buffer

    integer           :: i, j
    integer           :: type
    character(len=12) :: str1, str2
    
    if ( ascii ) then
      write(str1(1:12),'(i12)') nelems_final
      write(str2(1:12),'(i12)') (nodes_per_elem+1)*nelems_final
      write(ifile,'(a)') 'CELLS '//str1//' '//str2
      do j = istart, iend
        write(ifile,'(99i10)') nodes_per_elem, elemnode(j,:)-1
      end do
    else
      write(str1(1:12),'(i12)') nelems_final
      write(str2(1:12),'(i12)') (nodes_per_elem+1)*nelems_final
      buffer = lf//'CELLS '//str1//' '//str2//lf  ; write(ifile) trim(buffer)
      do j = istart, iend
        write(ifile) int(nodes_per_elem,4)
        write(ifile) (int(elemnode(j,i)-1,4),i=1,nodes_per_elem)
      end do
    end if

    if (nodes_per_elem==6) then
      type=13
    elseif (nodes_per_elem==8) then
      type=12
    else
      write(*,*) 'not implemented'
      stop
    end if

    if ( ascii ) then
      write(str1(1:12),'(i12)') nelems_final
      write(ifile,'(a)') 'CELL_TYPES'//str1
      write(ifile,'(4i10)') (type, i=1,nelems_final)
    else
      write(str1(1:12),'(i12)') nelems_final
      buffer = lf//'CELL_TYPES '//str1//lf; write(ifile) trim(buffer)
      write(ifile) (int(type,4), i=1,nelems_final)
    end if

  end subroutine write_cell


  
  subroutine write_header_scalar(ifile, nelems)
    !> header for scalar quantities (needed once for multiple quantities)
    integer, intent(in)             :: nelems, ifile
    character(len=80) :: buffer
    
    character(len=12) :: str1
    if (ascii) then
      write(ifile,'(a,i12)')'CELL_DATA ',nelems
    else
      write(str1, '(i12)') nelems
      buffer = lf//'CELL_DATA '//str1//lf; write(ifile) trim(buffer)
    end if
  end subroutine write_header_scalar
  
  subroutine write_scalar(name, ifile, istart, iend, values)
    !> values of scalars
    real*8 , intent(in)             :: values(:)
    integer, intent(in)             :: istart, iend, ifile
    character(len=*), intent(in)    :: name
    
    integer :: i

    if ( ascii ) then
      write(ifile,'(a)')'Scalars '//trim(name)//' float 1'
      write(ifile,'(a)')'LOOKUP_TABLE default'
      write(ifile,'(4es28.9)') (values(i), i=istart, iend)
    else
      buffer = lf//'Scalars '//trim(name)//' float 1'//lf;write(ifile) trim(buffer)
      buffer = lf//'LOOKUP_TABLE default'//lf;write(ifile) trim(buffer)
      write(ifile) (real(values(i),4), i=istart, iend)
    endif

  end subroutine write_scalar
  
  subroutine write_vector(name, ifile, istart, iend, v1, v2, v3)
    !> Write vector information with header
    real*8 , intent(in)             :: v1(:),v2(:),v3(:)
    integer, intent(in)             :: istart, iend, ifile
    character(len=*), intent(in)    :: name
    
    integer :: i
    if ( ascii ) then
      write(ifile,'(a,i12)')'VECTORS '//trim(name)//' float'
      write(ifile,'(4es28.9)') ((/v1(i),v2(i),v3(i)/), i=istart, iend)
    else
      buffer = lf//'VECTORS '//trim(name)//' float'//lf; write(ifile) trim(buffer)
      write(ifile) ( (/real(v1(i),4), real(v2(i),4), real(v3(i),4)/), i=istart, iend)
    endif

  end subroutine write_vector


  subroutine det_comp(elem_component, comp1, comp2, nelems, istart, iend, nelems_final)
    !> determine start and end index for different components
    integer, intent(in)             :: nelems, comp1, comp2
    integer, intent(in)             :: elem_component(:)
    integer, intent(out)            :: istart, iend, nelems_final

    character(len=80)               :: buffer
    integer           :: i
    istart = -1
    if (comp1 ==0) then
      nelems_final = nelems
      istart = 1; iend = nelems
    else
      nelems_final=0
      do i=1, nelems
        if (elem_component(i)>=comp1 .and. elem_component(i)<=comp2 ) then
          if (istart == -1) istart = i
          iend = i
          nelems_final = nelems_final+1
        end if
      end do
    end if
  end subroutine det_comp



subroutine write_header_pov(ifile)
    !> head for pov, viewing angles are hardcoded
  integer,              intent(in) :: ifile
  
  980 format(a)
  write(ifile,980) '#include "colors.inc"'
  write(ifile,980) '#include "textures.inc"'
  write(ifile,980) ''
  write(ifile,980) 'global_settings{'
  write(ifile,980) '    assumed_gamma   1.0 '
  write(ifile,980) '    max_trace_level 50'
  write(ifile,980) '}'
  write(ifile,980) ''
  write(ifile,980) 'camera{'
  write(ifile,980) '    angle 42'
  write(ifile,980) '    location<0.3,-11,4>'
  write(ifile,980) '    look_at <0.3,0,0>'
  write(ifile,980) '    up<0,0,6>'
  write(ifile,980) '    right<8,0,0>'
  write(ifile,980) '}'
  write(ifile,980) ''
  write(ifile,980) 'light_source{'
  write(ifile,980) '    <+0.5,-3,3>'
  write(ifile,980) '    color White'
  write(ifile,980) '    shadowless'
  write(ifile,980) '}'
  write(ifile,980) ''
  write(ifile,980) 'light_source{'
  write(ifile,980) '    <-0.3,+1.2,-3.0>'
  write(ifile,980) '    color rgb<0.5,0.5,0.5>'
  write(ifile,980) '    shadowless'
  write(ifile,980) '}'
  write(ifile,980) ''
  write(ifile,980) 'background{'
  write(ifile,980) '    White'
  write(ifile,980) '}'
  write(ifile,980) '  '
  
end subroutine write_header_pov



subroutine write_currdens_pov(ifile)
  !> write coordinates and scalar 
  integer,              intent(in) :: ifile
  
  integer :: i, j
  real*8, dimension(8) :: x, y, z
  real*8, dimension(4)  :: rgbt
  
  do i = 1, n_elems
    if ( .not. include_struct(elem_component(i)) ) cycle
    
    select case (povray_colormap)
    case (1)
      rgbt = colormap1(wall_curr_el(i), povray_wallcurr_min, povray_wallcurr_max)
    case default
      rgbt = colormap0(wall_curr_el(i), povray_wallcurr_min, povray_wallcurr_max)
    end select
    
    do j = 1, nodes_per_elem
      
      x(j) = xyznode(elemnode(i,j),1)
      y(j) = xyznode(elemnode(i,j),2)
      z(j) = xyznode(elemnode(i,j),3)
      
    end do
    phi = atan2(sum(y)/real(nodes_per_elem),sum(x)/real(nodes_per_elem))
    if (phi<0) phi = phi + 2* pi
    if (phi > phimax) cycle
    if (phi < phimin) cycle
    call write_hexahedron_as_triangles_pov(ifile, x, y, z, rgbt)
    
  end do

end subroutine write_currdens_pov



subroutine write_triangle_pov(ifile, x, y, z, rgbt)
  integer,              intent(in) :: ifile
  real*8, dimension(3), intent(in) :: x, y, z
  real*8, dimension(4), intent(in) :: rgbt
  
  990 format(a)
  991 format('  <',f15.6,',',f15.6,',',f15.6,'>, <',f15.6,',',f15.6,',',f15.6,'>, <',f15.6, &
    ',',f15.6,',',f15.6,'>')
  992 format('  pigment{color rgbt<',f15.6,','f15.6,','f15.6,','f15.6,'>}')
  write(ifile,990) 'triangle {'
  write(ifile,991) x(1), y(1), z(1), x(2), y(2), z(2), x(3), y(3), z(3)
  write(ifile,992) rgbt
  write(ifile,990) '}'

end subroutine write_triangle_pov



subroutine write_hexahedron_as_triangles_pov(ifile, x, y, z, rgbt)
  !> subdivide hexagon information into triangles
  integer,              intent(in) :: ifile
  real*8, dimension(8), intent(in) :: x, y, z
  real*8, dimension(4), intent(in) :: rgbt
  
  call write_triangle_pov(ifile, (/x(1), x(2), x(3)/), (/y(1), y(2), y(3)/), (/z(1), z(2), z(3)/), rgbt)
  call write_triangle_pov(ifile, (/x(1), x(3), x(4)/), (/y(1), y(3), y(4)/), (/z(1), z(3), z(4)/), rgbt)
  
  call write_triangle_pov(ifile, (/x(5), x(6), x(7)/), (/y(5), y(6), y(7)/), (/z(5), z(6), z(7)/), rgbt)
  call write_triangle_pov(ifile, (/x(5), x(7), x(8)/), (/y(5), y(7), y(8)/), (/z(5), z(7), z(8)/), rgbt)
  
  call write_triangle_pov(ifile, (/x(1), x(2), x(5)/), (/y(1), y(2), y(5)/), (/z(1), z(2), z(5)/), rgbt)
  call write_triangle_pov(ifile, (/x(2), x(5), x(6)/), (/y(2), y(5), y(6)/), (/z(2), z(5), z(6)/), rgbt)
  
  call write_triangle_pov(ifile, (/x(2), x(3), x(6)/), (/y(2), y(3), y(6)/), (/z(2), z(3), z(6)/), rgbt)
  call write_triangle_pov(ifile, (/x(3), x(6), x(7)/), (/y(3), y(6), y(7)/), (/z(3), z(6), z(7)/), rgbt)
  
  call write_triangle_pov(ifile, (/x(3), x(4), x(7)/), (/y(3), y(4), y(7)/), (/z(3), z(4), z(7)/), rgbt)
  call write_triangle_pov(ifile, (/x(4), x(7), x(8)/), (/y(4), y(7), y(8)/), (/z(4), z(7), z(8)/), rgbt)
  
  call write_triangle_pov(ifile, (/x(4), x(1), x(8)/), (/y(4), y(1), y(8)/), (/z(4), z(1), z(8)/), rgbt)
  call write_triangle_pov(ifile, (/x(1), x(8), x(5)/), (/y(1), y(8), y(5)/), (/z(1), z(8), z(5)/), rgbt)
  
end subroutine write_hexahedron_as_triangles_pov



function colormap0(v, vmin, vmax)
  real*8, dimension(4) :: colormap0
  real*8,               intent(in)  :: v, vmin, vmax
  
  colormap0 = (/ 1.d0, 1.d0, 1.d0, 0.0d0 /)
  
end function colormap0



function colormap1(v, vmin, vmax)
  real*8, dimension(4) :: colormap1
  real*8,               intent(in)  :: v, vmin, vmax
  
  real*8 :: v_norm
  
  if ( vmax<=vmin ) then
    colormap1 = (/ 0.d0, 0.d0, 0.d0, 1.d0 /) ! write invisible triangle
  else
    v_norm = min(max(v,vmin),vmax) / (vmax-vmin) ! crop to min/max range and normalize
    
    colormap1 = (/ 1.d0, 1.d0-v_norm, 1.d0-v_norm, 0.0d0 /)
  end if
  
end function colormap1


  
end program convert


