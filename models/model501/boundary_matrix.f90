subroutine boundary_matrix(vertex, element, nodes, xpoint2, xcase2, psi_axis, psi_bnd, Z_xpoint, ELM, RHS)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundaries of one element
! implements the natural boundary conditions
!---------------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(2)        ! the two nodes containing the boundary nodes

real*8     :: x_g(n_gauss), x_s(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss)
real*8     :: eq_ss(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

!real*8, dimension (:,:), pointer  :: ELM
!real*8, dimension (:)  , pointer  :: RHS
real*8     :: ELM(n_vertex_max*n_var*n_degrees*n_tor,n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: RHS(n_vertex_max*n_var*n_degrees*n_tor)

integer    :: vertex(2), i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: ws, xjac,  BigR, phi, eps_cyl
real*8     :: psi_axis, psi_bnd, Z_xpoint(2)
real*8     :: rhs_ij_6
real*8     :: psi_norm, theta, zeta, gamma_sheeth

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, Vpar0, r0, T0  
real*8     :: psi, psi_s, vpar, rho,  T   
real*8     :: amat_61, amat_65, amat_66, amat_67
logical    :: xpoint2

!theta = 0.5d0; zeta = 0.d0          ! Crank-Nicholson parameter
!theta = 1.0d0  ; zeta = 0.0d0       ! Euler scheme 
!theta = 1.0d0   ; zeta = 0.5d0      ! BDF2 (Gears) scheme

theta=time_evol_theta
zeta=time_evol_zeta

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0;  x_ss  = 0.d0; 
y_g  = 0.d0; y_s  = 0.d0;  y_ss  = 0.d0; 
eq_g = 0.d0; eq_s = 0.d0;  eq_ss = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0

do i=1,2
  
  do j=1,2

    do ms=1, n_gauss

      x_g(ms)  = x_g(ms)  + nodes(i)%x(1,j,1) * element%size(i,j) * H1(i,j,ms)
      x_s(ms)  = x_s(ms)  + nodes(i)%x(1,j,1) * element%size(i,j) * H1_s(i,j,ms)

      y_g(ms)  = y_g(ms)  + nodes(i)%x(1,j,2) * element%size(i,j) * H1(i,j,ms)
      y_s(ms)  = y_s(ms)  + nodes(i)%x(1,j,2) * element%size(i,j) * H1_s(i,j,ms)

      do mp=1,n_plane

        do k=1,n_var

          do in=1,n_tor

            eq_g(mp,k,ms)  = eq_g(mp,k,ms)  + nodes(i)%values(in,j,k) * element%size(i,j) * H1(i,j,ms)   * HZ(in,mp)

            eq_s(mp,k,ms)  = eq_s(mp,k,ms)  + nodes(i)%values(in,j,k) * element%size(i,j) * H1_s(i,j,ms) * HZ(in,mp)
	    
            eq_p(mp,k,ms)  = eq_p(mp,k,ms)  + nodes(i)%values(in,j,k) * element%size(i,j) * H1(i,j,ms)   * HZ_p(in,mp)

            eq_ss(mp,k,ms) = eq_ss(mp,k,ms) + nodes(i)%values(in,j,k) * element%size(i,j) * H1_ss(i,j,ms)* HZ(in,mp)

            delta_g(mp,k,ms) = delta_g(mp,k,ms) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H1(i,j,ms)   * HZ(in,mp)
            delta_s(mp,k,ms) = delta_s(mp,k,ms) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H1_s(i,j,ms) * HZ(in,mp)

          enddo
        enddo
      enddo

    enddo
  enddo
enddo


!gamma_sheeth = -3.d0

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

   ws = wgauss(ms)

   do mp = 1, n_plane

     ps0   = eq_g(mp,1,ms)
     ps0_s = eq_s(mp,1,ms)

     r0    = eq_g(mp,5,ms)
     T0    = eq_g(mp,6,ms)
     Vpar0 = eq_g(mp,7,ms)

     psi_norm = (ps0 - psi_axis)/(psi_bnd - psi_axis)
     if (xpoint2) then
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms) .lt. Z_xpoint(1)) .and. (xcase2 .ne. UPPER_XPOINT)) then
         psi_norm = 2.d0 - psi_norm
       endif
       if ((psi_norm .lt. 1.d0) .and. (y_g(ms) .gt. Z_xpoint(2)) .and. (xcase2 .ne. LOWER_XPOINT)) then
         psi_norm = 2.d0 - psi_norm
       endif
     endif


     do i=1,2                ! loop over nodes
     
       do j=1,2              ! loop over basis functions

         do im=1,n_tor

           index_ij = n_tor*n_var*n_degrees*(vertex(i)-1) + n_tor * n_var * (j-1) + im   ! index in the ELM matrix

           v   =  H1(i,j,ms) * element%size(vertex(i),j) * HZ(im,mp)

           rhs_ij_6 = v * gamma_sheath * r0 * T0 * vpar0 * ps0_s * tstep 

           ij6 = index_ij + 5*n_tor

           RHS(ij6) = RHS(ij6) + rhs_ij_6 * ws

           do k=1,2          ! loop over nodes

             do l=1,2        ! loop over basis functions

               do in = 1, n_tor

                 psi   = H1(k,l,ms)   * element%size(vertex(k),l) * HZ(in,mp)

                 psi_s = H1_s(k,l,ms) * element%size(vertex(k),l) * HZ(in,mp)

                 rho   = psi    ;    T   = psi   ;    vpar   = psi

                 index_kl = n_tor*n_var*n_degrees*(vertex(k)-1) + n_tor * n_var * (l-1) + in   ! index in the ELM matrix


                 amat_61 = - v * gamma_sheath * r0  * T0 * vpar0 * psi_s * theta * tstep 

                 amat_65 = - v * gamma_sheath * rho * T0 * vpar0 * ps0_s * theta * tstep 

                 amat_66 = - v * gamma_sheath * r0  * T  * vpar0 * ps0_s * theta * tstep 

                 amat_67 = - v * gamma_sheath * r0  * T0 * vpar  * ps0_s * theta * tstep 
		 

                 kl1 = index_kl
                 kl5 = index_kl + 4*n_tor
                 kl6 = index_kl + 5*n_tor
                 kl7 = index_kl + 6*n_tor

                 ELM(ij6,kl1) =  ELM(ij6,kl1) + ws * amat_61
                 ELM(ij6,kl5) =  ELM(ij6,kl5) + ws * amat_65
                 ELM(ij6,kl6) =  ELM(ij6,kl6) + ws * amat_66
                 ELM(ij6,kl7) =  ELM(ij6,kl7) + ws * amat_67

               enddo
             enddo
           enddo

         enddo
       enddo
     enddo
     
   enddo
enddo

return
end
