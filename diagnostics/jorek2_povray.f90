!> 3D visualizations of JOREK data with the raytracing program Povray.
!!
!! 2023 by Matthias Hoelzl and Nina Schwarz
!!
program jorek2_povray

use constants
use data_structure
use nodes_elements
use phys_module
use mod_import_restart
use mod_interp
use mod_basisfunctions
use mod_boundary
use equil_info

implicit none

integer, parameter :: ifile   = 53 !< File handle
logical, parameter :: surface = .true.
real*8,  parameter :: phimin  = -180.d0 / 360.d0 * 2.d0 * PI
real*8,  parameter :: phimax  =  +90.d0 / 360.d0 * 2.d0 * PI
integer, parameter :: n_phi   = 64
integer, parameter :: nsub    = 4
integer, parameter :: ivar    = 3
real*8,  parameter :: valmin  = -8.d0
real*8,  parameter :: valmax  = 8.d0

integer :: ierr, my_id, i_elm, i_s, i_t, iv, idof, itor, iplane, idir
real*8  :: s, t, phi, v1, v2
real*8  :: HH(nsub+1,nsub+1,4,n_degrees), H1(nsub+1,2,2), HZ(n_tor, n_phi)
real*8  :: val(nsub+1,nsub+1), vv(nsub+1, n_phi)
real*8  :: R(nsub+1,nsub+1),   RR(nsub+1, n_phi)
real*8  :: Z(nsub+1,nsub+1),   ZZ(nsub+1, n_phi)
real*8  :: x(nsub+1,nsub+1),   xx(nsub+1, n_phi)
real*8  :: y(nsub+1,nsub+1),   yy(nsub+1, n_phi)
real*8  :: rgbt1(4), rgbt2(4)
type(type_element)     :: element
type(type_bnd_element) :: bnd_element
type(type_node)        :: nodes(4)

! --- Initialize
my_id = 0
call initialise_parameters(my_id, "__NO_FILENAME__")
call det_modes()
call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)

open(ifile, file='jorek_3d.pov', form='formatted', status='replace')
call write_header_pov(ifile)

do i_s = 1, nsub+1
  s = real(i_s-1)/real(nsub)
  call basisfunctions1(s, H1(i_s, :, :))
end do

do i_s = 1, nsub+1
  s = real(i_s-1)/real(nsub)
  do i_t = 1, nsub+1
    t = real(i_t-1)/real(nsub)
    call basisfunctions(s, t, HH(i_s, i_t, :, :))
  end do
end do

do itor = 1, n_tor
  do iplane = 1, n_phi
    phi = phimin + real(iplane-1)/real(n_phi-1) * (phimax - phimin)
    if ( (itor == 1) .or. (mod(itor,2)==1) ) then
      HZ(itor,iplane) = cos(mode(itor)*phi)
    else
      HZ(itor,iplane) = sin(mode(itor)*phi)
    end if
  end do
end do

! --- Write the first and last poloidal cut
do iplane = 1, n_phi, n_phi-1
  phi = phimin + real(iplane-1)/real(n_phi-1) * (phimax - phimin)
  do i_elm = 1, element_list%n_elements
    element  = element_list%element(i_elm)
    nodes(:) = node_list%node(element%vertex(:))
    val(:,:) = 0.d0
    R(:,:)   = 0.d0
    Z(:,:)   = 0.d0
    do itor = 1, n_tor
      do i_s = 1, nsub+1
        do i_t = 1, nsub+1
          do iv = 1, 4
            do idof = 1, 4
              if ( itor == 1 ) then
                R  (i_s,i_t) = R  (i_s,i_t) + nodes(iv)%x(1,idof,1) * HH(i_s, i_t, iv, idof) * element%size(iv, idof)
                Z  (i_s,i_t) = Z  (i_s,i_t) + nodes(iv)%x(1,idof,2) * HH(i_s, i_t, iv, idof) * element%size(iv, idof)
              end if
              val(i_s,i_t) = val(i_s,i_t) + nodes(iv)%values(itor,idof,ivar) * HH(i_s, i_t, iv, idof) * element%size(iv, idof) * HZ(itor,iplane)
            end do
          end do
        end do
      end do
    end do
    
    x(:,:) =   R(:,:) * cos(phi)
    y(:,:) = - R(:,:) * sin(phi)
    
    do i_s = 1, nsub
      do i_t = 1, nsub
        v1 = ( val(i_s,i_t) + val(i_s+1,i_t) + val(i_s+1,i_t+1) ) / 3.d0
        v2 = ( val(i_s,i_t) + val(i_s,i_t+1) + val(i_s+1,i_t+1) ) / 3.d0
        rgbt1 = colormap1(v1, valmin, valmax)
        rgbt2 = colormap1(v2, valmin, valmax)
        call write_triangle_pov(ifile, &
          (/ x(i_s,i_t), x(i_s+1,i_t), x(i_s+1,i_t+1) /), &
          (/ y(i_s,i_t), y(i_s+1,i_t), y(i_s+1,i_t+1) /), &
          (/ z(i_s,i_t), z(i_s+1,i_t), z(i_s+1,i_t+1) /), &
          rgbt1)
        call write_triangle_pov(ifile, &
          (/ x(i_s,i_t), x(i_s,i_t+1), x(i_s+1,i_t+1) /), &
          (/ y(i_s,i_t), y(i_s,i_t+1), y(i_s+1,i_t+1) /), &
          (/ z(i_s,i_t), z(i_s,i_t+1), z(i_s+1,i_t+1) /), &
          rgbt2)
      end do
    end do
    
  end do
end do

! --- Write the boundary
do i_elm = 1, bnd_elm_list%n_bnd_elements
  bnd_element = bnd_elm_list%bnd_element(i_elm)
  nodes(1:2)  = node_list%node(bnd_element%vertex(:))
  
  RR(:,:) = 0.d0
  ZZ(:,:) = 0.d0
  vv(:,:) = 0.d0
  do i_s = 1, nsub+1
    do iplane = 1, n_phi
      phi = phimin + real(iplane-1)/real(n_phi-1) * (phimax - phimin)
      do iv = 1, 2
        do idof = 1, 2
          idir = bnd_element%direction(iv,idof)
          RR(i_s,iplane) = RR(i_s,iplane) + nodes(iv)%x(1,idir,1) * H1(i_s, iv, idof) * bnd_element%size(iv, idof)
          ZZ(i_s,iplane) = ZZ(i_s,iplane) + nodes(iv)%x(1,idir,2) * H1(i_s, iv, idof) * bnd_element%size(iv, idof)
          xx(i_s,iplane) =   RR(i_s,iplane) * cos(phi)
          yy(i_s,iplane) = - RR(i_s,iplane) * sin(phi)
          do itor = 1, n_tor
            vv(i_s,iplane) = vv(i_s,iplane) + nodes(iv)%values(itor,idir,ivar) * H1(i_s, iv, idof) * bnd_element%size(iv, idof) * HZ(itor,iplane)
          end do
        end do
      end do
    end do
  end do
  
  do i_s = 1, nsub
    do iplane = 1, n_phi-1
      v1 = ( vv(i_s,iplane) + vv(i_s+1,iplane) + vv(i_s+1,iplane+1) ) / 3.d0
      v2 = ( vv(i_s,iplane) + vv(i_s,iplane+1) + vv(i_s+1,iplane+1) ) / 3.d0
      rgbt1 = colormap1(v1, valmin, valmax)
      rgbt2 = colormap1(v2, valmin, valmax)
      call write_triangle_pov(ifile, &
        (/ xx(i_s,iplane), xx(i_s+1,iplane), xx(i_s+1,iplane+1) /), &
        (/ yy(i_s,iplane), yy(i_s+1,iplane), yy(i_s+1,iplane+1) /), &
        (/ zz(i_s,iplane), zz(i_s+1,iplane), zz(i_s+1,iplane+1) /), &
        rgbt1)
      call write_triangle_pov(ifile, &
        (/ xx(i_s,iplane), xx(i_s,iplane+1), xx(i_s+1,iplane+1) /), &
        (/ yy(i_s,iplane), yy(i_s,iplane+1), yy(i_s+1,iplane+1) /), &
        (/ zz(i_s,iplane), zz(i_s,iplane+1), zz(i_s+1,iplane+1) /), &
        rgbt2)
    end do
  end do

end do

close(ifile)



contains



subroutine write_header_pov(ifile)
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
  write(ifile,980) '    location<0.3,-11,-4>'
  write(ifile,980) '    look_at <0.3,0,0>'
  write(ifile,980) '    up<0,0,-6>'
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



function colormap1(v, vmin, vmax)
  real*8, dimension(4) :: colormap1
  real*8,               intent(in)  :: v, vmin, vmax
  
  real*8 :: v_norm
  
  if ( vmax<=vmin ) then
    colormap1 = (/ 0.d0, 0.d0, 0.d0, 1.d0 /) ! write invisible triangle
  else
    v_norm = min(max(v,vmin),vmax) / (vmax-vmin) ! crop to min/max range and normalize
    
    colormap1 = (/ min(1.d0,2.d0-2.d0*v_norm), 1.d0-2.d0*abs(v_norm-0.5d0), min(1.d0,2.d0*v_norm), 0.0d0 /)
  end if
  
end function colormap1


  
end program jorek2_povray
