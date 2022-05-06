subroutine element_matrix_projection(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
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
use mod_interp
use nodes_elements
use mod_sources

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss), x_s(n_gauss,n_gauss), x_t(n_gauss,n_gauss)
real*8     :: y_g(n_gauss,n_gauss), y_s(n_gauss,n_gauss), y_t(n_gauss,n_gauss)
real*8     :: factor(n_gauss,n_gauss)
real*8     :: eq_g(n_gauss,n_gauss),   eq_s(n_gauss,n_gauss),   eq_t(n_gauss,n_gauss)
real*8     :: eq2_g(n_gauss,n_gauss),  eq2_s(n_gauss,n_gauss),  eq2_t(n_gauss,n_gauss)
real*8     :: ELM(n_vertex_max*n_degrees,n_vertex_max*n_degrees), RHS(n_vertex_max*n_degrees)

real*8     :: xjac, wst
real*8     :: v, psi, rhs_ij
integer    :: ms, mt, i, j, k, l, index_ij, index_kl, itype, ivar_in, ivar_out, i_harm
real*8     :: F_prof        ,dF_dpsi      ,dF_dz      , dF_dpsi2      ,dF_dz2       ,dF_dpsi_dz
real*8     :: zFFprime      ,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi2,dFFprime_dz2 ,dFFprime_dpsi_dz
real*8     :: zn,  dn_dpsi,  dn_dpsi2,  dn_dz,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi2_dz,  dn_dpsi_dz2
real*8     :: zT,  dT_dpsi,  dT_dpsi2,  dT_dz,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi2_dz,  dT_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: zrn, drn_dpsi, drn_dpsi2, drn_dz, drn_dz2, drn_dpsi_dz, drn_dpsi3, drn_dpsi2_dz, drn_dpsi_dz2
real*8     :: zV,  dV_dpsi,  dV_dpsi2,  dV_dz,  dV_dz2,  dV_dpsi_dz,  dV_dpsi3,  dV_dpsi2_dz,  dV_dpsi_dz2
real*8     :: var_RHS
real*8     :: R_out,Z_out,s_out,t_out
integer    :: ielm_out,ifail
real*8     :: dd1,dd2,dd3,dd4,dd5

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

       if (ivar_out .eq. 710) then
#ifdef fullmhd
         eq2_g(ms,mt)  = eq2_g(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H(i,j,ms,mt)
         eq2_s(ms,mt)  = eq2_s(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H_s(i,j,ms,mt)
         eq2_t(ms,mt)  = eq2_t(ms,mt)  + nodes(i)%Fprof_eq(j) * element%size(i,j) * H_t(i,j,ms,mt)
#endif
       else
         eq2_g(ms,mt)  = eq2_g(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H(i,j,ms,mt)
         eq2_s(ms,mt)  = eq2_s(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H_s(i,j,ms,mt)
         eq2_t(ms,mt)  = eq2_t(ms,mt)  + nodes(i)%values(i_harm,j,ivar_out) * element%size(i,j) * H_t(i,j,ms,mt)
       endif

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
    
    if (ivar_out .eq. 710) then
      ! --- note: no need to use psi_axis_init, psi_bnd_init etc., since this routine should only be called at t=0
      call F_profile      (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           F_prof        ,dF_dpsi      ,dF_dz      , dF_dpsi2      ,dF_dz2       ,dF_dpsi_dz , &
                           zFFprime      ,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi2,dFFprime_dz2 ,dFFprime_dpsi_dz)
      var_RHS = F_prof
    endif
    if ( (ivar_out .eq. var_rho) .or. (ivar_out .eq. var_zj) ) then
      call density        (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
      if (ivar_out .eq. var_rho) var_RHS = zn
    endif
    if ( (ivar_out .eq. var_T) .or. ( (ivar_out .eq. var_zj) .and. (.not. with_TiTe)) ) then
      call temperature    (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
      if (ivar_out .eq. var_T) var_RHS = zT
    endif
    if ( (ivar_out .eq. var_Ti) .or. ( (ivar_out .eq. var_zj) .and. with_TiTe) ) then
      call temperature_i  (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
      if (ivar_out .eq. var_Ti) var_RHS = zTi
    endif
    if ( (ivar_out .eq. var_Te) .or. ( (ivar_out .eq. var_zj) .and. with_TiTe) ) then
      call temperature_e  (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
      if (ivar_out .eq. var_Te) var_RHS = zTe
    endif
    if (ivar_out .eq. var_zj) then
      call FFprime        (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz, .true.)
    endif
    if (ivar_out .eq. var_Vpar) then
      call velocity       (xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zV,dV_dpsi,dV_dz,dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
      var_RHS = zV
    endif
    if (ivar_out .eq. var_rhon) then
      call neutral_density(xpoint, xcase, y_g(ms,mt), ES%Z_xpoint, eq_g(ms,mt), ES%psi_axis, ES%psi_bnd, &
                           zrn,drn_dpsi,drn_dz,drn_dpsi2,drn_dz2,drn_dpsi_dz,drn_dpsi3,drn_dpsi_dz2, drn_dpsi2_dz)
      var_RHS = zrn
    endif
    
    ! --- Special case for psi, we interpolate from the old grid copied into grid_xpoint_data module
    if (ivar_out .eq. var_psi) then
      call find_RZ(node_list,element_list,x_g(ms,mt),y_g(ms,mt),R_out,Z_out,ielm_out,s_out,t_out,ifail)
      if (ifail .ne. 0) then ! if we can't find the point on the previous grid, use whatever was set before (should be the eqdsk value) 
        !write(*,'(A,2f18.7)')'Projection Warning: location not found on old grid',x_g(ms,mt),y_g(ms,mt)
        var_RHS = eq2_g(ms,mt)
      else
        call interp(node_list,element_list,ielm_out,1,1,s_out,t_out,psi,dd1,dd2,dd3,dd4,dd5)
        var_RHS = psi
      endif
    endif

    ! --- Special case for current
    if (ivar_out .eq. var_zj) then
      if (with_TiTe) then
        zT = zTi + zTe
        dT_dpsi = dTi_dpsi + dTe_dpsi
      endif
      var_RHS = zFFprime - x_g(ms,mt)**2 * (dn_dpsi * zT + zn * dT_dpsi)
    endif
 
    do i=1,n_vertex_max
 
      do j=1,n_degrees
 
        index_ij = (i-1)*n_degrees + j
 
        v   = h(i,j,ms,mt)  * element%size(i,j)
 
        rhs_ij = + var_RHS
 
        RHS(index_ij) = RHS(index_ij) - v * rhs_ij       * factor(ms,mt) * xjac * wst
        RHS(index_ij) = RHS(index_ij) + v * eq2_g(ms,mt) * factor(ms,mt) * xjac * wst    ! solve for perturbation only
 
        do k=1,n_vertex_max
 
          do l=1,n_degrees
 
            psi   = h(k,l,ms,mt)  * element%size(k,l)
 
            index_kl = (k-1)*n_degrees + l
 
            ELM(index_ij,index_kl) =  ELM(index_ij,index_kl) - psi * v * factor(ms,mt) * xjac * wst
 
          enddo
        enddo
 
      enddo
    enddo
 
  enddo
enddo

return
end
