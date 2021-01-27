subroutine plot_solution(node_list,element_list,ivar,iharm,iangle,label)
!-----------------------------------------------------------------------
! plots the 2D solution
!-----------------------------------------------------------------------
use constants
use tr_module 
use mod_parameters
use data_structure
use basis_at_gaussian
use phys_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list

real*8             :: psi_max, psi_min, R_min, R_max, Z_min, Z_max, zc(3)
real*8,allocatable :: xp(:,:)
real*8             :: x_g(4,4), x_s(4,4), x_t(4,4), y_g(4,4), y_s(4,4), y_t(4,4)
real*8             :: PS(4,4),  PSX(4,4), PSY(4,4), xjac, error, error2, error3
real*8             :: x_hel, y_hel, ps_solov, psx_solov, psy_solov, psxy_solov
integer            :: i, kv, kf, ms, mt, nc, iv, ivar, iharm, nlab, in, iangle
character*(*)      :: label
logical            :: all_harmonics

if ( .not. write_ps ) then
  write(*,*) 'Jorek postscript deactivated. Skip plot_solutions'
  return
endif
write(*,*) '****************************************************'
write(*,'(A,A,2i5)') ' * plotting solution : ',label,iharm,iangle
write(*,*) '****************************************************'
write(*,*) ' elements : ',element_list%n_elements
write(*,*) ' nodes    : ',node_list%n_nodes

all_harmonics = .true.
if ( (iharm .ge. 0) .and. (iharm .le. n_tor) ) all_harmonics = .false.

if (all_harmonics) write(*,*) ' total ',ivar,iangle

psi_max = -1.e20
psi_min = +1.e20

do i=1,element_list%n_elements
 if(element_list%element(i)%sons(1).eq.0) then
 x_g = 0. ; x_s = 0.; x_t = 0.
 y_g = 0. ; y_s = 0.; y_t = 0.
 PS = 0.  ; PSX = 0.; PSY = 0.

 do ms = 1, 4          ! 4 Gaussian points
   do mt = 1, 4        ! 4 Gaussian points

     do kf = 1, 4       ! 4 basis functions
       do kv = 1, 4     ! 4 vertices

          iv = element_list%element(i)%vertex(kv)

          x_g(ms,mt) = x_g(ms,mt) + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

       enddo
     enddo

   enddo
 enddo

 do ms = 1, 4          ! 4 Gaussian points
   do mt = 1, 4        ! 4 Gaussian points

     xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

     do kf = 1, 4       ! 4 basis functions
       do kv = 1, 4     ! 4 vertices

         iv = element_list%element(i)%vertex(kv)

         if ( .not. all_harmonics ) then

           PS(ms,mt)  = PS(ms,mt)  + node_list%node(iv)%values(iharm,kf,ivar) &
                        * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

         else

           PS(ms,mt)  = PS(ms,mt)  + node_list%node(iv)%values(1,kf,ivar)     &
                        * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt) * HZ(1,iangle)

           do in=2, n_tor

             PS(ms,mt)  = PS(ms,mt) + node_list%node(iv)%values(in,kf,ivar) &
                        * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt) * HZ(in,iangle)
           enddo

         endif
       enddo
     enddo

   enddo
 enddo

 psi_max = max(psi_max,maxval(PS))
 psi_min = min(psi_min,minval(PS))
endif
enddo

if ((abs(psi_max) .lt. 1.e-30) .and. (abs(psi_min) .lt. 1.e-30)) then
  write(*,*) ' empty plot',psi_min,psi_max
  return
endif

write(*,*) ' min/max : ',psi_min,psi_max

call tr_allocate(xp,1,node_list%n_nodes,1,n_dim,"xp",CAT_GRID)

do i=1,node_list%n_nodes
 xp(i,1:n_dim) = node_list%node(i)%x(1,1,1:n_dim)
enddo
R_max = maxval(xp(:,1))
R_min = minval(xp(:,1))
Z_max = maxval(xp(:,2))
Z_min = minval(xp(:,2))

!R_min  = 1.1 * R_min - 0.1 * R_max
!Z_min  = 1.1 * Z_min - 0.1 * Z_max

nlab = len(label)
call nframe(21,11,1,R_min,R_max,Z_min,Z_max,LABEL,nlab,'R [m]',5,'Z [m]',5)

!call lplot(21,11,791,xp(:,1),xp(:,2),-node_list%n_nodes,1,'Nodes',5,'X',1,'Y',1)

nc = 3
zc(1) = -max(abs(psi_max),abs(psi_min))
zc(2) = 0.
zc(3) = abs(zc(1))

error = 0.

do i=1,element_list%n_elements
 if(element_list%element(i)%sons(1).eq.0) then
 x_g = 0. ; x_s = 0.; x_t = 0.
 y_g = 0. ; y_s = 0.; y_t = 0.
 PS = 0.  ; PSX = 0.; PSY = 0.

 do ms = 1, 4          ! 4 Gaussian points
   do mt = 1, 4        ! 4 Gaussian points

     do kf = 1, 4       ! 4 basis functions
       do kv = 1, 4     ! 4 vertices

          iv = element_list%element(i)%vertex(kv)

          x_g(ms,mt) = x_g(ms,mt) + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)

       enddo
     enddo

   enddo
 enddo

 do ms = 1, 4          ! 4 Gaussian points
   do mt = 1, 4        ! 4 Gaussian points

     xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

     do kf = 1, 4       ! 4 basis functions
       do kv = 1, 4     ! 4 vertices

         iv = element_list%element(i)%vertex(kv)

         if ( .not. all_harmonics ) then

           PS(ms,mt)  = PS(ms,mt)  + node_list%node(iv)%values(iharm,kf,ivar) &
                        * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)

         else

           PS(ms,mt)  = PS(ms,mt)  + node_list%node(iv)%values(1,kf,ivar)     &
                        * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt) * HZ(1,iangle)

           do in=2, n_tor

             PS(ms,mt)  = PS(ms,mt) + node_list%node(iv)%values(in,kf,ivar) &
                        * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt) * HZ(in,iangle)
           enddo

         endif
       enddo
     enddo

   enddo
 enddo

 call cplotm(1,1,-2,x_g,y_g,4,-4,1,1,PS,4,zc,nc,'Solution',8,'R [m]',5,'Z [m]',5)
 endif
enddo

call cplotm(1,1,-1,x_g,y_g,4,-4,1,1,PS,4,zc,nc,'Solution',8,'R [m]',5,'Z [m]',5)
call lincol(0)

if (allocated(xp)) call tr_deallocate(xp,"xp",CAT_GRID)

return
end

