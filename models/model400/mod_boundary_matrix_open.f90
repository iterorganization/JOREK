module mod_boundary_matrix_open
  implicit none
contains

subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max)
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
use grid_xpoint_data

implicit none

type (type_element)          :: element
type (type_node)             :: nodes(2)        ! the two nodes containing the boundary nodes
type (type_strategic_points) :: stpts
integer,     intent(in)      :: i_tor_min   
integer,     intent(in)      :: i_tor_max   

real*8     :: x_g(n_gauss), x_s(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss), eq_ss(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

real*8     :: ELM(n_vertex_max*n_var*(n_order+1)*n_tor,n_vertex_max*n_var*(n_order+1)*n_tor)
real*8     :: RHS(n_vertex_max*n_var*(n_order+1)*n_tor)

integer    :: vertex(2), direction(2), i, j, j2, ms, mt, mp, k, l, l2, index_ij, index_kl, index, xcase2
integer    :: in, im, n_tor_local
integer    :: ij1, ij2, ij3, ij4, ij5, ij6, ij7, ij8
integer    :: kl1, kl2, kl3, kl4, kl5, kl6, kl7, kl8
real*8     :: ws, xjac,  BigR, Z, phi, eps_cyl
real*8     :: psi_axis, R_axis, Z_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
real*8     :: rhs_ij_5, rhs_ij_6, rhs_ij_7, rhs_ij_8
real*8     :: theta, zeta

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, Vpar0, r0, Ti0, Te0  
real*8     :: psi, psi_s, vpar, rho,  Ti,  Te
real*8     :: alpha, normal, DL, Dwall, gas_puff
real*8     :: amat_51, amat_55, amat_57
real*8     :: amat_61, amat_65, amat_66, amat_67
real*8     :: amat_81, amat_85, amat_88, amat_87
real*8     :: element_size_ij, element_size_kl
logical    :: xpoint2

!theta = 0.5d0; zeta = 0.d0          ! Crank-Nicholson parameter
!theta = 1.0d0  ; zeta = 0.0d0       ! Euler scheme 
!theta = 1.0d0   ; zeta = 0.5d0       ! BDF2 (Gears) scheme

theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0;  x_ss  = 0.d0; 
y_g  = 0.d0; y_s  = 0.d0;  y_ss  = 0.d0; 
eq_g = 0.d0; eq_s = 0.d0;  eq_ss = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0;

Dwall    = 0.d0!0.3d0
gas_puff = 0.d-8

do i=1,2

  do j=1,2

    j2 = direction(j)

    element_size_ij = element%size(vertex(i),j2)

    do ms=1, n_gauss

      x_g(ms)  = x_g(ms)  + nodes(i)%x(1,j2,1) * element_size_ij * H1(i,j,ms)
      x_s(ms)  = x_s(ms)  + nodes(i)%x(1,j2,1) * element_size_ij * H1_s(i,j,ms)

      y_g(ms)  = y_g(ms)  + nodes(i)%x(1,j2,2) * element_size_ij * H1(i,j,ms)
      y_s(ms)  = y_s(ms)  + nodes(i)%x(1,j2,2) * element_size_ij * H1_s(i,j,ms)

      do mp=1,n_plane

        do k=1,n_var

          do in=1,n_tor

            eq_g(mp,k,ms)  = eq_g(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp)
            eq_s(mp,k,ms)  = eq_s(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ(in,mp)
            eq_p(mp,k,ms)  = eq_p(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ_p(in,mp)
            eq_ss(mp,k,ms) = eq_ss(mp,k,ms) + nodes(i)%values(in,j2,k) * element_size_ij * H1_ss(i,j,ms)* HZ(in,mp)

            delta_g(mp,k,ms) = delta_g(mp,k,ms) + nodes(i)%deltas(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp)
            delta_s(mp,k,ms) = delta_s(mp,k,ms) + nodes(i)%deltas(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ(in,mp)

          enddo
        enddo
      enddo

    enddo
  enddo
enddo

! changes deltas for variable time steps
delta_g = delta_g * tstep / tstep_prev
delta_s = delta_s * tstep / tstep_prev

n_tor_local = i_tor_max - i_tor_min + 1
!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

   ws = wgauss(ms)

   DL   = sqrt(x_s(ms)**2 + y_s(ms)**2)
   BigR = x_g(ms)
   Z    = y_g(ms)

   do mp = 1, n_plane

     ps0   = eq_g(mp,1,ms)
     ps0_s = eq_s(mp,1,ms)             ! why not absolute value for normal orientation?

     r0    = eq_g(mp,5,ms)
     Ti0   = eq_g(mp,6,ms)
     Te0   = eq_g(mp,8,ms)
     Vpar0 = eq_g(mp,7,ms)

     if (xcase2 .eq. 1) then
       if (BigR .lt. (stpts%RRightCorn_LowerInnerLeg+stpts%RLeftCorn_LowerOuterLeg)/2.0) then
         normal = -1.0
       else
         normal = +1.0
       endif
     endif
     if (xcase2 .eq. 2) then
       if (BigR .lt. (stpts%RRightCorn_UpperInnerLeg+stpts%RLeftCorn_UpperOuterLeg)/2.0) then
         normal = -1.0
       else
         normal = +1.0
       endif
     endif
     if (xcase2 .eq. 3) then
       if (Z .gt. (Z_xpoint(1) + Z_xpoint(2))/2.0) then
         if (BigR .lt. (stpts%RRightCorn_UpperInnerLeg+stpts%RLeftCorn_UpperOuterLeg)/2.0) then
           normal = -1.0
         else
           normal = +1.0
         endif
       else
         if (BigR .lt. (stpts%RRightCorn_LowerInnerLeg+stpts%RLeftCorn_LowerOuterLeg)/2.0) then
           normal = -1.0
         else
           normal = +1.0
         endif
       endif
     endif

     do i=1,2                ! loop over nodes

       do j=1,2              ! loop over basis functions

         j2 = direction(j)
         element_size_ij = element%size(vertex(i),j2)

         do im= i_tor_min, i_tor_max

           index_ij = n_tor_local*n_var*(n_order+1)*(vertex(i)-1) + n_tor_local * n_var * (j2-1) + im - i_tor_min + 1  ! index in the ELM matrix

           v   =  H1(i,j,ms) * element_size_ij * HZ(im,mp)         ! test function

!diffuse perp flux only (unless we change the integration by parts of the density equation)
           rhs_ij_5 = + v * density_reflection * r0 * vpar0 * ps0_s * normal * tstep  &      ! right hand side equation 5
                      + v * gas_puff * DL * BigR                             * tstep  &
                      - v * Dwall * r0 * r0 * DL                             * tstep !* BigR

           rhs_ij_6 = - v * (gamma_sheath -1.d0) * r0 * Ti0 * vpar0 * ps0_s * normal * tstep  ! right hand side equation 6
           rhs_ij_8 = - v * (gamma_sheath -1.d0) * r0 * Te0 * vpar0 * ps0_s * normal * tstep  ! right hand side equation 6

           ij5 = index_ij + 4*n_tor_local                                          ! local index in element matrix
           ij6 = index_ij + 5*n_tor_local                                          ! local index in element matrix
           ij8 = index_ij + 7*n_tor_local                                          ! local index in element matrix

           RHS(ij5) = RHS(ij5) + rhs_ij_5 * ws                               ! add to element RHS
           RHS(ij6) = RHS(ij6) + rhs_ij_6 * ws                               ! add to element RHS
           RHS(ij8) = RHS(ij8) + rhs_ij_8 * ws                               ! add to element RHS

           do k=1,2                                                          ! loop over nodes

             do l=1,2                                                        ! loop over basis functions

               l2 = direction(l)
               element_size_kl = element%size(vertex(k),l2)

               do in =  i_tor_min, i_tor_max                                              ! loop over toroidal harmonics

                 psi   = H1(k,l,ms)   * element_size_kl * HZ(in,mp)

                 psi_s = H1_s(k,l,ms) * element_size_kl * HZ(in,mp)

                 rho   = psi    ;    Ti   = psi   ;    Te   = psi   ;    vpar   = psi

!diffuse perp flux only (unless we change the integration by parts of the density equation)

                 amat_51 = - v * density_reflection * r0  * vpar0 * psi_s * normal * theta * tstep 
                 amat_55 = - v * density_reflection * rho * vpar0 * ps0_s * normal * theta * tstep &
                           + v * Dwall * 2.d0 * r0 * rho * DL                      * theta * tstep ! * BigR
                 amat_57 = - v * density_reflection * r0  * vpar  * ps0_s * normal * theta * tstep 

                 amat_61 = + v * (gamma_sheath-1.d0) * r0  * Ti0 * vpar0 * psi_s * normal * theta * tstep 
                 amat_65 = + v * (gamma_sheath-1.d0) * rho * Ti0 * vpar0 * ps0_s * normal * theta * tstep &
                           + v * (gamma_sheath-1.d0) * rho * Ti0 * 0.00  * DL    * normal * theta * tstep 
                 amat_66 = + v * (gamma_sheath-1.d0) * r0  * Ti  * vpar0 * ps0_s * normal * theta * tstep &
                           + v * (gamma_sheath-1.d0) * r0  * Ti  * 0.00  * DL    * normal * theta * tstep 
                 amat_67 = + v * (gamma_sheath-1.d0) * r0  * Ti0 * vpar  * ps0_s * normal * theta * tstep 

                 amat_81 = + v * (gamma_sheath-1.d0) * r0  * Te0 * vpar0 * psi_s * normal * theta * tstep 
                 amat_85 = + v * (gamma_sheath-1.d0) * rho * Te0 * vpar0 * ps0_s * normal * theta * tstep &
                           + v * (gamma_sheath-1.d0) * rho * Te0 * 0.00  * DL    * normal * theta * tstep 
                 amat_88 = + v * (gamma_sheath-1.d0) * r0  * Te  * vpar0 * ps0_s * normal * theta * tstep &
                           + v * (gamma_sheath-1.d0) * r0  * Te  * 0.00  * DL    * normal * theta * tstep 
                 amat_87 = + v * (gamma_sheath-1.d0) * r0  * Te0 * vpar  * ps0_s * normal * theta * tstep 

                 index_kl = n_tor_local*n_var*(n_order+1)*(vertex(k)-1) + n_tor_local * n_var * (l2-1) + in - i_tor_min + 1  ! index in the ELM matrix

                 kl1 = index_kl
                 kl5 = index_kl + 4*n_tor_local
                 kl6 = index_kl + 5*n_tor_local
                 kl7 = index_kl + 6*n_tor_local
                 kl8 = index_kl + 7*n_tor_local

                 ELM(ij5,kl1) =  ELM(ij5,kl1) + ws * amat_51
                 ELM(ij5,kl5) =  ELM(ij5,kl5) + ws * amat_55
                 ELM(ij5,kl7) =  ELM(ij5,kl7) + ws * amat_57

                 ELM(ij6,kl1) =  ELM(ij6,kl1) + ws * amat_61
                 ELM(ij6,kl5) =  ELM(ij6,kl5) + ws * amat_65
                 ELM(ij6,kl6) =  ELM(ij6,kl6) + ws * amat_66
                 ELM(ij6,kl7) =  ELM(ij6,kl7) + ws * amat_67

                 ELM(ij8,kl1) =  ELM(ij8,kl1) + ws * amat_81
                 ELM(ij8,kl5) =  ELM(ij8,kl5) + ws * amat_85
                 ELM(ij8,kl8) =  ELM(ij8,kl8) + ws * amat_88
                 ELM(ij8,kl7) =  ELM(ij8,kl7) + ws * amat_87

               enddo
             enddo
           enddo

         enddo
       enddo
     enddo

   enddo
enddo

return
end subroutine

end module mod_boundary_matrix_open
