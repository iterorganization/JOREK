subroutine export_POV(node_list,element_list,ivar,iharm)
!**************************************************************************
! write an input file for POVRAY using the bezier patch object            *
!**************************************************************************
use mod_parameters
use data_structure
use basis_at_gaussian
implicit none

type(type_node_list)    :: node_list
type(type_element_list) :: element_list

real*8 :: x11(3),x14(3),x44(3),x41(3)
real*8 :: u11(3),u14(3),u44(3),u41(3)
real*8 :: v11(3),v14(3),v44(3),v41(3)
real*8 :: w11(3),w14(3),w44(3),w41(3)

real*8  :: x_min, x_max, y_min, y_max, z_min, z_max, x_centre, y_centre, scale, rr, gg, bb
integer :: ivar, iharm, iv1, iv2, iv3, iv4, i
character*11 :: filename

x_min = 1.d10; x_max = -1.d10; y_min = 1.d10; y_max = -1.d10;  z_min = 1.d10; z_max = -1.d10
do i=1,node_list%n_nodes
 x_min = min(x_min,node_list%node(i)%x(1,1,1))
 x_max = max(x_max,node_list%node(i)%x(1,1,1))
 y_min = min(y_min,node_list%node(i)%x(1,1,2))
 y_max = max(y_max,node_list%node(i)%x(1,1,2))
 z_min = min(z_min,real(node_list%node(i)%values(iharm,1,ivar)))
 z_max = max(z_max,real(node_list%node(i)%values(iharm,1,ivar)))
enddo

x_centre = (x_min + x_max) /2.d0
y_centre = (y_min + y_max) /2.d0


write(*,*) ' export to POVRAY : ',element_list%n_elements
write(*,*) ' x_min, x_max, y_min, y_max : ',x_min, x_max, y_min, y_max
write(*,*) ' z_min, z_max               : ',z_min, z_max

write(filename,'(A6,i1,A4)') 'jorek_',ivar,'.pov'

open(42,file=filename)

write(42,*) '#include "shapes.inc"   '
write(42,*) '#include "colors.inc"   '
write(42,*) '#include "textures.inc" '
write(42,*) '#include "glass.inc" '

write(42,*) ' plane {  z, -100 hollow on   texture { pigment { color red 0.4 green 0.4 blue 0.4 } } } '
write(42,'(A,f8.3,A,f8.3,A,f8.3,A1,f8.3,A)') ' camera { location  <',x_centre,',',y_centre, &
                                            ', 4.0>   look_at <',x_centre,',',y_centre,', 0.0>} '
write(42,*) ' light_source { <-2, 2, 5> colour White } '
write(42,*) ' light_source { <2, 2, 5> colour White } '
write(42,*) ' light_source { <2, -2, 5> colour White } '
write(42,*) ' light_source { <-2, -2, 5> colour White } '

scale = 2.d0 / (z_max - z_min)

do i=1,element_list%n_elements

 write(42,*) 'bicubic_patch { type 1 flatness 0.  u_steps 3  v_steps 3'

! rr = (1.d0+mod(i,3)) / 3.d0
! gg = (1.d0+mod(i,5)) / 5.d0
! bb = (1.d0+mod(i,7)) / 7.d0

 rr = 0.5
 gg = 0.5
 bb = 0.5

 iv1 = element_list%element(i)%vertex(1)
 iv2 = element_list%element(i)%vertex(2)
 iv3 = element_list%element(i)%vertex(3)
 iv4 = element_list%element(i)%vertex(4)

 x11 = 0.d0; x14 = 0.d0; x44 = 0.d0; x41 = 0.d0
 u11 = 0.d0; u14 = 0.d0; u44 = 0.d0; u41 = 0.d0
 v11 = 0.d0; v14 = 0.d0; v44 = 0.d0; v41 = 0.d0
 w11 = 0.d0; w14 = 0.d0; w44 = 0.d0; w41 = 0.d0

 x11(1:2) = node_list%node(iv1)%x(1,1,1:2)
 x14(1:2) = node_list%node(iv2)%x(1,1,1:2)
 x44(1:2) = node_list%node(iv3)%x(1,1,1:2)
 x41(1:2) = node_list%node(iv4)%x(1,1,1:2)

 do iharm=1,n_tor
   x11(3) = x11(3) + scale * node_list%node(iv1)%values(iharm,1,ivar) * HZ(iharm,1)
   x14(3) = x14(3) + scale * node_list%node(iv2)%values(iharm,1,ivar) * HZ(iharm,1)
   x44(3) = x44(3) + scale * node_list%node(iv3)%values(iharm,1,ivar) * HZ(iharm,1)
   x41(3) = x41(3) + scale * node_list%node(iv4)%values(iharm,1,ivar) * HZ(iharm,1)
 enddo

 u11(1:2) = node_list%node(iv1)%x(1,2,1:2) * element_list%element(i)%size(1,2)
 u14(1:2) = node_list%node(iv2)%x(1,2,1:2) * element_list%element(i)%size(2,2)
 u44(1:2) = node_list%node(iv3)%x(1,2,1:2) * element_list%element(i)%size(3,2)
 u41(1:2) = node_list%node(iv4)%x(1,2,1:2) * element_list%element(i)%size(4,2)

 do iharm=1,n_tor
   u11(3) = u11(3) + scale * node_list%node(iv1)%values(iharm,2,ivar) * element_list%element(i)%size(1,2) * HZ(iharm,1)
   u14(3) = u14(3) + scale * node_list%node(iv2)%values(iharm,2,ivar) * element_list%element(i)%size(2,2) * HZ(iharm,1)
   u44(3) = u44(3) + scale * node_list%node(iv3)%values(iharm,2,ivar) * element_list%element(i)%size(3,2) * HZ(iharm,1)
   u41(3) = u41(3) + scale * node_list%node(iv4)%values(iharm,2,ivar) * element_list%element(i)%size(4,2) * HZ(iharm,1)
 enddo

 v11(1:2) = node_list%node(iv1)%x(1,3,1:2) * element_list%element(i)%size(1,3)
 v14(1:2) = node_list%node(iv2)%x(1,3,1:2) * element_list%element(i)%size(2,3)
 v44(1:2) = node_list%node(iv3)%x(1,3,1:2) * element_list%element(i)%size(3,3)
 v41(1:2) = node_list%node(iv4)%x(1,3,1:2) * element_list%element(i)%size(4,3)

 do iharm=1,n_tor
   v11(3) = v11(3) + scale * node_list%node(iv1)%values(iharm,3,ivar) * element_list%element(i)%size(1,3) * HZ(iharm,1)
   v14(3) = v14(3) + scale * node_list%node(iv2)%values(iharm,3,ivar) * element_list%element(i)%size(2,3) * HZ(iharm,1)
   v44(3) = v44(3) + scale * node_list%node(iv3)%values(iharm,3,ivar) * element_list%element(i)%size(3,3) * HZ(iharm,1)
   v41(3) = v41(3) + scale * node_list%node(iv4)%values(iharm,3,ivar) * element_list%element(i)%size(4,3) * HZ(iharm,1)
 enddo

 w11(1:2) = node_list%node(iv1)%x(1,4,1:2) * element_list%element(i)%size(1,4)
 w14(1:2) = node_list%node(iv2)%x(1,4,1:2) * element_list%element(i)%size(2,4)
 w44(1:2) = node_list%node(iv3)%x(1,4,1:2) * element_list%element(i)%size(3,4)
 w41(1:2) = node_list%node(iv4)%x(1,4,1:2) * element_list%element(i)%size(4,4)

 do iharm=1,n_tor
   w11(3) = w11(3) + scale * node_list%node(iv1)%values(iharm,4,ivar) * element_list%element(i)%size(1,4) * HZ(iharm,1)
   w14(3) = w14(3) + scale * node_list%node(iv2)%values(iharm,4,ivar) * element_list%element(i)%size(2,4) * HZ(iharm,1)
   w44(3) = w44(3) + scale * node_list%node(iv3)%values(iharm,4,ivar) * element_list%element(i)%size(3,4) * HZ(iharm,1)
   w41(3) = w41(3) + scale * node_list%node(iv4)%values(iharm,4,ivar) * element_list%element(i)%size(4,4) * HZ(iharm,1)
 enddo

 w11 = x11 + u11 + v11 + w11
 w14 = x14 + u14 + v14 + w14
 w44 = x44 + u44 + v44 + w44
 w41 = x41 + u41 + v41 + w41

 u11 = x11 + u11
 u14 = x14 + u14
 u44 = x44 + u44
 u41 = x41 + u41

 v11 = x11 + v11
 v14 = x14 + v14
 v44 = x44 + v44
 v41 = x41 + v41


 write(42,'(12(A,e12.4),A)')   '<',x11(1),',',x11(2),',',x11(3),'>, <',u11(1),',',u11(2),',',u11(3), &
                          '>,<',u14(1),',',u14(2),',',u14(3),'>, <',x14(1),',',x14(2),',',x14(3),'>'

 write(42,'(12(A,e12.4),A)')   '<',v11(1),',',v11(2),',',v11(3),'>, <',w11(1),',',w11(2),',',w11(3), &
                          '>,<',w14(1),',',w14(2),',',w14(3),'>, <',v14(1),',',v14(2),',',v14(3),'>'

 write(42,'(12(A,e12.4),A)')   '<',v41(1),',',v41(2),',',v41(3),'>, <',w41(1),',',w41(2),',',w41(3), &
                          '>,<',w44(1),',',w44(2),',',w44(3),'>, <',v44(1),',',v44(2),',',v44(3),'>'

 write(42,'(12(A,e12.4),A)')   '<',x41(1),',',x41(2),',',x41(3),'>, <',u41(1),',',u41(2),',',u41(3), &
                          '>,<',u44(1),',',u44(2),',',u44(3),'>, <',x44(1),',',x44(2),',',x44(3),'>'

 write(42,'(A,f5.1,A,f5.1,A,f5.1,A)') 'texture { pigment { color red',rr,' green ',gg,' blue ',bb, &
                                    '} finish { roughness 0.001 specular 0.6 diffuse 0.4} } }'

!endif
enddo


close(42)

return
end