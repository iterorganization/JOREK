subroutine element_matrix_710_equi(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use equil_info, only: ES
use phys_module, only: xpoint, xcase
use mod_F_profile

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: factor(n_gauss,n_gauss)
real*8     :: eq_g(n_gauss,n_gauss),   eq_s(n_gauss,n_gauss),   eq_t(n_gauss,n_gauss)
real*8     :: eq2_g(n_gauss,n_gauss),  eq2_s(n_gauss,n_gauss),  eq2_t(n_gauss,n_gauss)
real*8     :: ELM(n_vertex_max*(n_order+1),n_vertex_max*(n_order+1)), RHS(n_vertex_max*(n_order+1))

real*8     :: xjac, wst
real*8     :: v, psi, rhs_ij
integer    :: ms, mt, i, j, k, l, index_ij, index_kl, itype, ivar_in, ivar_out, i_harm
real*8     :: F_prof        ,dF_dpsi      ,dF_dz      , dF_dpsi2      ,dF_dz2       ,dF_dpsi_dz
real*8     :: zFFprime      ,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi2,dFFprime_dz2 ,dFFprime_dpsi_dz
#ifdef fullmhd

ELM=0.d0
RHS=0.d0

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g(:,:)   = 0.d0; x_s(:,:)   = 0.d0; x_t(:,:)   = 0.d0;
y_g(:,:)   = 0.d0; y_s(:,:)   = 0.d0; y_t(:,:)   = 0.d0;
eq_g(:,:)  = 0.d0; eq_s(:,:)  = 0.d0; eq_t(:,:)  = 0.d0;
eq2_g(:,:) = 0.d0; eq2_s(:,:) = 0.d0; eq2_t(:,:) = 0.d0;

do i=1,n_vertex_max
 do j=1,n_order+1
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

       eq2_g(ms,mt)  = eq2_g(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H(i,j,ms,mt)
       eq2_s(ms,mt)  = eq2_s(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H_s(i,j,ms,mt)
       eq2_t(ms,mt)  = eq2_t(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H_t(i,j,ms,mt)

     enddo
   enddo
 enddo
enddo


factor =  x_g                                   ! Poisson

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
   
   ! --- note: no need to use psi_axis_init, psi_bnd_init etc., since this routine should only be called at t=0
   call F_profile   (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                     F_prof        ,dF_dpsi      ,dF_dz      , dF_dpsi2      ,dF_dz2       ,dF_dpsi_dz , &
                     zFFprime      ,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi2,dFFprime_dz2 ,dFFprime_dpsi_dz)

   do i=1,n_vertex_max

     do j=1,n_order+1

       index_ij = (i-1)*(n_order+1) + j

       v   = h(i,j,ms,mt)  * element%size(i,j)

       rhs_ij = + F_prof

       RHS(index_ij) = RHS(index_ij) - v * rhs_ij       * factor(ms,mt) * xjac * wst
       RHS(index_ij) = RHS(index_ij) + v * eq2_g(ms,mt) * factor(ms,mt) * xjac * wst    ! solve for perturbation only

       do k=1,n_vertex_max

         do l=1,n_order+1

           psi   = h(k,l,ms,mt)  * element%size(k,l)

           index_kl = (k-1)*(n_order+1) + l

           ELM(index_ij,index_kl) =  ELM(index_ij,index_kl) - psi * v * factor(ms,mt) * xjac * wst

         enddo
       enddo

     enddo
   enddo

 enddo
enddo

#endif
return
end
