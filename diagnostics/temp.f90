subroutine temp(node_list,element_list,A_tem,A_den,A_jen,A_jec,A_jec1,A_jec2)
!---------------------------------------------------------------
! subroutine to drop global (broken down by n_tor) 
! temperature, density, current, and eccd current using
! live data routines in communication folder
!---------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use mod_parameters

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8     :: eq_g(n_var,n_gauss,n_gauss), eq_s(n_var,n_gauss,n_gauss), eq_t(n_var,n_gauss,n_gauss)

integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements
real*8     :: A_tem(n_tor), A_den(n_tor), xjac, BigR, wst
real*8     :: A_jen(n_tor), A_jec(n_tor), A_jec1(n_tor),A_jec2(n_tor)
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t
real*8     :: t0, t0_x, t0_y, t0_p, t0_s, t0_t
real*8     :: p0, p0_s, p0_t, p0_x, p0_y, zj0
real*8     :: jec, jec1, jec2

A_tem = 0.d0
A_den = 0.d0
A_jen = 0.d0
A_jec = 0.d0
A_jec1 = 0.d0
A_jec2 = 0.d0


do ife =1,  element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:) = 0.d0;    x_s(:,:) = 0.d0;    x_t(:,:) = 0.d0;
  y_g(:,:) = 0.d0;    y_s(:,:) = 0.;      y_t(:,:) = 0.d0;
  eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_degrees
      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

        enddo
      enddo
    enddo
  enddo

  do in=1,n_tor

    eq_g(:,:,:) = 0.d0; eq_s(:,:,:) = 0.d0; eq_t(:,:,:) = 0.d0;

    do i=1,n_vertex_max
      do j=1,n_degrees
        do ms=1, n_gauss
          do mt=1, n_gauss

            do k=1,n_var
              eq_g(k,ms,mt)  = eq_g(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)
              eq_s(k,ms,mt)  = eq_s(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
              eq_t(k,ms,mt)  = eq_t(k,ms,mt)  + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)
            enddo
        
          enddo
        enddo
      enddo
    enddo



!--------------------------------------------------- sum over the Gaussian integration points
    do ms=1, n_gauss
      do mt=1, n_gauss
        wst = wgauss(ms)*wgauss(mt)
        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

!density definition
     r0    = abs(eq_g(5,ms,mt))

! temperature definition
     t0    = abs(eq_g(6,ms,mt))

! pressure definition
     p0    = r0 * t0

! general current
     zj0   = eq_g(3,ms,mt)

! eccd current
   jec  = 0.d0
   if(jorek_model == 305) then
     jec=eq_g(n_var,ms,mt) ! the eccd current should always be the last equation
   elseif(jorek_model == 306) then
     jec1=eq_g(n_var-1,ms,mt)
     jec2=eq_g(n_var,ms,mt)
     jec =eq_g(n_var,ms,mt)+eq_g(n_var-1,ms,mt) ! the eccd current is the sum of two currents
   else
     jec = p0  ! when not applying eccd, this last array is global pressure
   endif
        A_tem(in) = A_tem(in) + t0* xjac * wst
        A_den(in) = A_den(in) + r0* xjac * wst
        A_jen(in) = A_jen(in) + zj0* xjac * wst
        if(jorek_model == 305) then
           A_jec(in) = A_jec(in) + jec* xjac * wst
        endif
        if(jorek_model == 306) then
           A_jec1(in) = A_jec1(in) + jec1* xjac * wst
           A_jec2(in) = A_jec2(in) + jec2* xjac * wst
           A_jec(in) = A_jec(in) + jec* xjac * wst
        endif
      enddo
    enddo
  enddo

enddo

return
end subroutine temp
