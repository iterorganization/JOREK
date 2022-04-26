subroutine element_matrix_GS(xpoint,xcase,Z_xpoint,psi_axis,psi_bnd,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: factor(n_gauss,n_gauss)
real*8     :: eq_g(n_gauss,n_gauss),  eq_s(n_gauss,n_gauss),  eq_t(n_gauss,n_gauss)
real*8     :: eq2_g(n_gauss,n_gauss), eq2_s(n_gauss,n_gauss), eq2_t(n_gauss,n_gauss)
real*8     :: ELM(n_vertex_max*n_degrees,n_vertex_max*n_degrees), RHS(n_vertex_max*n_degrees)

real*8     :: xjac, wst
real*8     :: ps0_x, ps0_y, v, v_x, v_y, psi, psi_x, psi_y, rhs_ij
integer    :: ms, mt, i, j, k, l, index_ij, index_kl, itype, ivar_in, ivar_out, i_harm,xcase
logical    :: xpoint
real*8     :: Z_xpoint(2),psi_axis,psi_bnd,dj_dpsi,dj_dz
real*8     :: zn,dn_dpsi,dn_dz, ddn_dpsi, ddn_dz, ddn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz
real*8     :: zT,dT_dpsi,dT_dz, ddT_dpsi, ddT_dz, ddT_dpsi_dz, dT_dpsi3, dT_dpsi_dz2, dT_dpsi2_dz
real*8     :: zTi,dTi_dpsi,dTi_dz, ddTi_dpsi, ddTi_dz, ddTi_dpsi_dz, dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz
real*8     :: zTe,dTe_dpsi,dTe_dz, ddTe_dpsi, ddTe_dz, ddTe_dpsi_dz, dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz
real*8     :: ddFFprime_dpsi_dz, zFFprime, dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi2,dFFprime_dz2


ELM=0.d0
RHS=0.d0

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g(:,:)   = 0.d0; x_s(:,:)   = 0.d0; x_t(:,:)   = 0.d0;
y_g(:,:)   = 0.d0; y_s(:,:)   = 0.d0; y_t(:,:)   = 0.d0;
eq_g(:,:)  = 0.d0; eq_s(:,:)  = 0.d0; eq_t(:,:)  = 0.d0;
eq2_g(:,:) = 0.d0; eq2_s(:,:) = 0.d0; eq2_t(:,:) = 0.d0;

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

       eq_g(ms,mt)  = eq_g(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * element%size(i,j) * H(i,j,ms,mt)
       eq_s(ms,mt)  = eq_s(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * element%size(i,j) * H_s(i,j,ms,mt)
       eq_t(ms,mt)  = eq_t(ms,mt)  + nodes(i)%values(i_harm,j,ivar_in) * element%size(i,j) * H_t(i,j,ms,mt)

       eq2_g(ms,mt)  = eq2_g(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H(i,j,ms,mt)
       eq2_s(ms,mt)  = eq2_s(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H_s(i,j,ms,mt)
       eq2_t(ms,mt)  = eq2_t(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H_t(i,j,ms,mt)

     enddo
   enddo
 enddo
enddo


factor = 1.d0/ x_g                    ! Grad-Shafranov

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

  do mt=1, n_gauss

    call density(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq2_g(ms,mt),psi_axis,psi_bnd, &
                zn,dn_dpsi,dn_dz,ddn_dpsi,ddn_dz,ddn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)

    if (with_TiTe) then  
       call temperature_i(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq2_g(ms,mt),psi_axis,psi_bnd, &
    			zTi,dTi_dpsi,dTi_dz,ddTi_dpsi,ddTi_dz,ddTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)
       call temperature_e(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq2_g(ms,mt),psi_axis,psi_bnd, &
    			zTe,dTe_dpsi,dTe_dz,ddTe_dpsi,ddTe_dz,ddTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)
       zT = zTi + zTe
       dT_dpsi = dTi_dpsi + dTe_dpsi	  
    else  
       call temperature(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq2_g(ms,mt),psi_axis,psi_bnd, &
    			zT,dT_dpsi,dT_dz,ddT_dpsi,ddT_dz,ddT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
    endif

    call FFprime(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq2_g(ms,mt),psi_axis,psi_bnd, &
                zFFprime, dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,ddFFprime_dpsi_dz, .true.)

    wst = wgauss(ms)*wgauss(mt)

    xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

    ps0_x = (   y_t(ms,mt) * eq2_s(ms,mt) - y_s(ms,mt) * eq2_t(ms,mt) ) / xjac
    ps0_y = ( - x_t(ms,mt) * eq2_s(ms,mt) + x_s(ms,mt) * eq2_t(ms,mt) ) / xjac

    do i=1,n_vertex_max

      do j=1,n_degrees

        index_ij = (i-1)*n_degrees + j

        v   = h(i,j,ms,mt)  * element%size(i,j)
        v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac
        v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac

        rhs_ij =  zFFprime / x_g(ms,mt) - (zn * dT_dpsi + dn_dpsi * zT) * x_g(ms,mt)

        RHS(index_ij) = RHS(index_ij) + v * rhs_ij  * xjac * wst

        do k=1,n_vertex_max

          do l=1,n_degrees

            psi   = h(k,l,ms,mt)  * element%size(k,l)
            psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) * element%size(k,l) / xjac
            psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) ) * element%size(k,l) / xjac

            index_kl = (k-1)*n_degrees + l

            ELM(index_ij,index_kl) =  ELM(index_ij,index_kl) - (psi_x * v_x + psi_y * v_y) * factor(ms,mt) * xjac * wst  !&
	   
!	    - 0.5d0 *( dFFprime_dpsi/x_g(ms,mt) - (zn*ddT_dpsi + ddn_dpsi*zT + 2.*dn_dpsi*dT_dpsi)*x_g(ms,mt)) *v*psi*xjac*wst

         enddo
       enddo

     enddo
   enddo

  enddo
enddo

return
end
