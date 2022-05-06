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
use corr_neg
use mod_interp
use diffusivities, only: get_dperp, get_zkperp

implicit none

integer,               intent(in)     :: i_tor_min   
integer,               intent(in)     :: i_tor_max   
! --- Routine variables
type (type_element),    intent(in)    :: element
type (type_node),       intent(inout) :: nodes(n_vertex_max)
integer,                intent(inout) :: vertex(2)
integer,                intent(in)    :: direction(2), xcase2
logical,                intent(in)    :: xpoint2
real*8,                 intent(in)    :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
real*8,                 intent(inout) :: ELM(n_vertex_max*n_var*n_degrees*n_tor,n_vertex_max*n_var*n_degrees*n_tor)
real*8,                 intent(inout) :: RHS(n_vertex_max*n_var*n_degrees*n_tor)

! --- Internal variables
real*8     :: x_g(n_gauss), x_s(n_gauss), x_t(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_t(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss), eq_ss(n_plane,n_var,n_gauss)
real*8     :: eq_t(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

integer    :: direction_perp(2), i, j, ms, mt, mp, k, l, index_ij, index_kl, index
integer    :: j2, l2, j3, l3, is, n_tor_local
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: ws, xjac, BigR, phi, eps_cyl, Btot, Btot_psi
real*8     :: rhs_ij_5, rhs_ij_6, rhs_ij_7
real*8     :: theta, zeta, Zbig, factor
real*8     :: R_mid, Z_mid, R_cnt, Z_cnt, normal_direction(2)

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, ps0_t, ps0_x, ps0_y, Vpar0, T0, T0_s, T0_t, T0_x, T0_y, T0_p, cs0, T_corr  
real*8     :: psi, psi_s, psi_t, psi_x, psi_y, vpar, T, cs_T, bdotn, bdotn2, B0_R, B0_Z, BB2
real*8     :: r0, r0_s, r0_t, r0_p, r0_x, r0_y, rho, rho_s, rho_t, rho_x, rho_y
real*8     :: amat_51, amat_55, amat_57
real*8     :: amat_61, amat_65, amat_66, amat_67
real*8     :: amat_71, amat_76, amat_77

real*8     :: alpha, R_inside, Z_inside, normal_sign, normal_sign2, normal_sign3, Dwall
real*8     :: element_size_ij, element_size_kl, element_size_perp, grad_t(2), normal(2)
real*8     :: threshold_cs, transition_cs

type (type_node)         :: tmp_node


theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_ss  = 0.d0; 
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_ss  = 0.d0; 
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_ss = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; 

direction_perp(1) = 6 / direction(2)     ! =3 if direction(2)=2, =2 if direction(2)=3
direction_perp(2) = 4

R_mid = sum(nodes(1:2)%x(1,1,1)) / 2.d0     ! mid point on boundary (approx.)
Z_mid = sum(nodes(1:2)%x(1,1,2)) / 2.d0
R_cnt = sum(nodes(1:4)%x(1,1,1)) / 4.d0     ! center point within element (approx.)
Z_cnt = sum(nodes(1:4)%x(1,1,2)) / 4.d0

normal_direction = (/R_mid - R_cnt, Z_mid - Z_cnt /) / norm2((/R_mid - R_cnt, Z_mid - Z_cnt /))


do i=1,2    ! sum over 2 verices
  
  do j=1,2  ! sum over two basis functions

    j2 = direction(j)
    element_size_ij = element%size(vertex(i),j2)
 
    j3 = direction_perp(j)
    element_size_perp = - element%size(vertex(i),direction_perp(1)) * 3.d0

    if ((vertex(1)*vertex(2) .eq. 2)) then
      element_size_perp = + element%size(vertex(i),direction_perp(1)) * 3.d0
    endif
    do ms=1, n_gauss

      x_g(ms)  = x_g(ms)  + nodes(i)%x(1,j2,1) * element_size_ij * H1(i,j,ms)
      x_s(ms)  = x_s(ms)  + nodes(i)%x(1,j2,1) * element_size_ij * H1_s(i,j,ms)
      x_t(ms)  = x_t(ms)  + nodes(i)%x(1,j3,1) * element_size_ij * H1(i,j,ms)   * element_size_perp

      y_g(ms)  = y_g(ms)  + nodes(i)%x(1,j2,2) * element_size_ij * H1(i,j,ms)
      y_s(ms)  = y_s(ms)  + nodes(i)%x(1,j2,2) * element_size_ij * H1_s(i,j,ms)
      y_t(ms)  = y_t(ms)  + nodes(i)%x(1,j3,2) * element_size_ij * H1(i,j,ms)   * element_size_perp

      do mp=1,n_plane

        do k=1,n_var

          do in=1,n_tor

            eq_g(mp,k,ms)  = eq_g(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp)

            eq_s(mp,k,ms)  = eq_s(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ(in,mp)

            eq_t(mp,k,ms)  = eq_t(mp,k,ms)  + nodes(i)%values(in,j3,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp) * element_size_perp

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

  xjac = x_s(ms)*y_t(ms) - x_t(ms)*y_s(ms)
  BigR = x_g(ms)

  grad_t = (/ - y_s(ms),   x_s(ms) /) / xjac

  normal = dot_product(grad_t,normal_direction) * grad_t      ! outward pointing normal
  normal = normal / norm2(normal)


  do mp = 1, n_plane

    ps0   = eq_g(mp,1,ms)
    ps0_s = eq_s(mp,1,ms)
    ps0_t = eq_t(mp,1,ms) 
    ps0_x = (   y_t(ms) * ps0_s - y_s(ms) * ps0_t ) / xjac
    ps0_y = ( - x_t(ms) * ps0_s + x_s(ms) * ps0_t ) / xjac

    r0   = eq_g(mp,5,ms)
    r0_s = eq_s(mp,5,ms)
    r0_t = eq_t(mp,5,ms)
    r0_p = eq_p(mp,5,ms)
    r0_x = (   y_t(ms) * r0_s - y_s(ms) * r0_t ) / xjac
    r0_y = ( - x_t(ms) * r0_s + x_s(ms) * r0_t ) / xjac

    T0     = eq_g(mp,6,ms)
    T0_s   = eq_s(mp,6,ms)
    T0_t   = eq_t(mp,6,ms)
    T0_p   = eq_p(mp,6,ms)
    T0_x   = (   y_t(ms) * T0_s - y_s(ms) * T0_t ) / xjac
    T0_y   = ( - x_t(ms) * T0_s + x_s(ms) * T0_t ) / xjac
    T_corr = max(T0,max(T_min,1.d-8))

    Vpar0 = eq_g(mp,7,ms)

    cs0   = sqrt(gamma*T_corr)

    B0_R =   ps0_y / x_g(ms)
    B0_Z = - ps0_x / x_g(ms)
    Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR
    BB2 = Btot**2
    !bdotn = (+ ps0_y * normal(1) - ps0_x * normal(2)) / x_g(ms) / Btot
    bdotn = (+ ps0_y * normal(1) - ps0_x * normal(2)) / sqrt(ps0_x**2 + ps0_y**2)
    !write(*,*)'check:',abs(bdotn)
    normal_sign  = bdotn / abs(bdotn)
    normal_sign3 = ps0_s / abs(ps0_s) * normal_sign
    
    ! --- transition factor between +cs0 and -cs0
    threshold_cs  = 5.d-1 !0.02d0
    transition_cs = 1.d-3 !0.016d0
    factor = 1.d0  * (0.5d0 + 0.5d0 * tanh((abs(bdotn) - threshold_cs)/transition_cs))**2
    Zbig   = 1.d20
    if ( (nodes(1)%boundary .eq. 9) .or. (nodes(2)%boundary .eq. 9) .or. (nodes(1)%boundary .eq. 3) .or. (nodes(2)%boundary .eq. 3) ) &
      Zbig   = 1.d20 * (0.5d0 + 0.5d0 * tanh((abs(bdotn) - threshold_cs)/transition_cs))**2

    do i=1,2                ! loop over nodes
    
      do j=1,2              ! loop over basis functions
    
        j2 = direction(j)
        element_size_ij = element%size(vertex(i),j2)

        do im=i_tor_min, i_tor_max

          index_ij = n_tor_local*n_var*n_degrees*(vertex(i)-1) + n_tor_local * n_var * (j2-1) + im - i_tor_min + 1  ! index in the ELM matrix

          v   =  H1(i,j,ms) * element_size_ij * HZ(im,mp)         ! test function

          rhs_ij_5 = + v * density_reflection * r0 * vpar0 * ps0_s * normal_sign3 * tstep          ! right hand side equation 5

          rhs_ij_6 = - v * (gamma_sheath -1.d0) * r0 * T0 * vpar0 * ps0_s * normal_sign3 * tstep   ! right hand side equation 6
          
          !rhs_ij_7 = - v * (vpar0 * Btot * normal_sign - cs0 * factor) * Zbig

          ij5 = index_ij + 4*n_tor_local                                          ! local index in element matrix
          ij6 = index_ij + 5*n_tor_local                                          ! local index in element matrix
          !ij7 = index_ij + 6*n_tor                                          ! local index in element matrix


          RHS(ij5) = RHS(ij5) + rhs_ij_5 * ws                               ! add to element RHS
          RHS(ij6) = RHS(ij6) + rhs_ij_6 * ws                               ! add to element RHS
          !RHS(ij7) = RHS(ij7) + rhs_ij_7 * ws                               ! add to element RHS
                   
          do k=1,2                                                          ! loop over nodes

            do l=1,2                                                        ! loop over basis functions
    
              l2 = direction(l)
              element_size_kl = element%size(vertex(k),l2)
 
              l3 = direction_perp(j)
              element_size_perp = - element%size(vertex(k),direction_perp(1)) * 3.d0

              do in = i_tor_min, i_tor_max                                              ! loop over toroidal harmonics

                psi   = H1(k,l,ms)   * element_size_kl * HZ(in,mp)
                psi_s = H1_s(k,l,ms) * element_size_kl * HZ(in,mp)
                psi_t = H1(k,l,ms)   * element_size_kl * HZ(in,mp) * element_size_perp
                psi_x = (   y_t(ms) * psi_s - y_s(ms) * psi_t ) / xjac
                psi_y = ( - x_t(ms) * psi_s + x_s(ms) * psi_t ) / xjac

                rho   = psi   ; T   = psi ; vpar   = psi
                rho_s = psi_s ;
                rho_t = psi_t ;
                rho_x = psi_x ;
                rho_y = psi_y ;

                cs_T  = gamma * T / (2.d0 * cs0)

                Btot_psi = (ps0_x*psi_x + ps0_y*psi_y) / sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

                amat_51 = - v * density_reflection * r0  * vpar0 * psi_s * normal_sign3 * theta * tstep 
                amat_55 = - v * density_reflection * rho * vpar0 * ps0_s * normal_sign3 * theta * tstep 
                amat_57 = - v * density_reflection * r0  * vpar  * ps0_s * normal_sign3 * theta * tstep 

                amat_61 = + v * (gamma_sheath-1.d0) * r0  * T0 * vpar0 * psi_s * normal_sign3 * theta * tstep 
                amat_65 = + v * (gamma_sheath-1.d0) * rho * T0 * vpar0 * ps0_s * normal_sign3 * theta * tstep 
                amat_66 = + v * (gamma_sheath-1.d0) * r0  * T  * vpar0 * ps0_s * normal_sign3 * theta * tstep 
                amat_67 = + v * (gamma_sheath-1.d0) * r0  * T0 * vpar  * ps0_s * normal_sign3 * theta * tstep 

                !amat_71 =   v * (vpar0 * Btot_psi * normal_sign) * Zbig
                !amat_76 =   v * ( - cs_T) * factor               * Zbig
                !amat_77 =   v * (vpar  * Btot     * normal_sign) * Zbig


                index_kl = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l2-1) + in - i_tor_min + 1  ! index in the ELM matrix
                
                kl1 = index_kl
                kl5 = index_kl + 4*n_tor_local
                kl6 = index_kl + 5*n_tor_local
                kl7 = index_kl + 6*n_tor_local

                ELM(ij5,kl1) =  ELM(ij5,kl1) + ws * amat_51
                ELM(ij5,kl5) =  ELM(ij5,kl5) + ws * amat_55
                ELM(ij5,kl7) =  ELM(ij5,kl7) + ws * amat_57

                ELM(ij6,kl1) =  ELM(ij6,kl1) + ws * amat_61
                ELM(ij6,kl5) =  ELM(ij6,kl5) + ws * amat_65
                ELM(ij6,kl6) =  ELM(ij6,kl6) + ws * amat_66
                ELM(ij6,kl7) =  ELM(ij6,kl7) + ws * amat_67

                !ELM(ij7,kl1) =  ELM(ij7,kl1) + ws * amat_71
                !ELM(ij7,kl6) =  ELM(ij7,kl6) + ws * amat_76
                !ELM(ij7,kl7) =  ELM(ij7,kl7) + ws * amat_77

              enddo
            enddo
          enddo

        enddo
      enddo
    enddo
    
  enddo
enddo

return
end subroutine boundary_matrix_open
end module mod_boundary_matrix_open


