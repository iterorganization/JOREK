subroutine integral_current(node_list,element_list,psi_axis, psi_bnd, xpoint2, xcase2, z_xpoint, current)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss)
real*8     :: ps_g(n_gauss,n_gauss)

integer    :: i, j, ms, mt, iv, inode, ife, n_elements, xcase2
real*8     :: zn,dn_dpsi,dn_dz,ddn_dpsi,ddn_dz,ddn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8     :: zT,dT_dpsi,dT_dz,ddT_dpsi,ddT_dz,ddT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz
real*8     :: zTi,zTe,dTi_dpsi,dTe_dpsi,dTi_dz,dTe_dz,ddTi_dpsi,ddTe_dpsi,ddTi_dz,ddTe_dz
real*8     :: ddTi_dpsi_dz,ddTe_dpsi_dz,dTi_dpsi3,dTe_dpsi3,dTi_dpsi_dz2,dTe_dpsi_dz2,dTi_dpsi2_dz,dTe_dpsi2_dz
real*8     :: zFFprime, dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,ddFFprime_dpsi_dz
real*8     :: current, xjac, BigR, Z_xpoint(2), psi_axis, psi_bnd, wst
logical    :: xpoint2

current = 0.d0

do ife =1,  element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)   = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0;
  y_g(:,:)   = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0;
  ps_g(:,:)  = 0.d0

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

          ps_g(ms,mt)  = ps_g(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H(i,j,ms,mt)

        enddo
      enddo
    enddo
  enddo

!--------------------------------------------------- sum over the Gaussian integration points

  do ms=1, n_gauss

    do mt=1, n_gauss

      call density(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, &
                   zn,dn_dpsi,dn_dz,ddn_dpsi,ddn_dz,ddn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)

      if (with_TiTe) then
        
        call temperature_i(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, &
                   zTi,dTi_dpsi,dTi_dz,ddTi_dpsi,ddTi_dz,ddTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)
        	       
        call temperature_e(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, &
                   zTe,dTe_dpsi,dTe_dz,ddTe_dpsi,ddTe_dz,ddTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)

        zT = zTi + zTe
        dT_dpsi = dTi_dpsi + dTe_dpsi
      
      else
        
        call temperature(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, &
                   zT,dT_dpsi,dT_dz,ddT_dpsi,ddT_dz,ddT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)

      endif

      call FFprime(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, ps_g(ms,mt),psi_axis,psi_bnd, &
                   zFFprime, dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,ddFFprime_dpsi_dz, .true.)
                   
      wst = wgauss(ms)*wgauss(mt)

      xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
      BigR = x_g(ms,mt)
           
      current = current + (zFFprime / x_g(ms,mt) - (zn * dT_dpsi + dn_dpsi * zT) * x_g(ms,mt)) * wst * xjac

    enddo
  enddo
enddo

current = current / (4.d-7 * PI)

write(*,'(A,f8.5,A)') ' current : ',current/1.e6,' MA'
return
end
