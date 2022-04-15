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

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)        ! the two nodes containing the boundary nodes
integer, intent(in)   :: i_tor_min   
integer, intent(in)   :: i_tor_max   

real*8     :: x_g(n_gauss), x_s(n_gauss), x_t(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_t(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss)
real*8     :: eq_t(n_plane,n_var,n_gauss), eq_ss(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

real*8     :: ELM(n_vertex_max*n_var*n_degrees*n_tor,n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: RHS(n_vertex_max*n_var*n_degrees*n_tor)

integer    :: vertex(2), direction(2), direction_perp(2)
integer    :: i, j, j2, j3, ms, mt, mp, k, l, l2, l3, index_ij, index_kl, index, xcase2, is
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, ij8, kl1, kl2, kl3, kl4, kl5, kl6, kl7, kl8
real*8     :: ws, xjac,  dl, BigR, phi, eps_cyl, Btot
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
real*8     :: rhs_ij_5, rhs_ij_6, rhs_ij_7, rhs_ij_8
real*8     :: theta, zeta, Zbig, BB2, bdotn, factor, psi_ss, vpar_ss
real*8     :: R_inside, Z_inside, R_mid, Z_mid, R_cnt, Z_cnt, normal(2), normal_direction(2)
real*8     :: normal_sign, normal_sign3

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, ps0_t, ps0_x, ps0_y, Vpar0, r0_corr, T0_corr, cs0  
real*8     :: psi, psi_s, psi_t, vpar, T, u0_s, u_s, cs_T
real*8     :: T0, T0_s, T0_t, T0_x, T0_y, T0_p
real*8     :: r0, r0_s, r0_t, r0_p, r0_x, r0_y, rho, rho_s, rho_t, rho_x, rho_y
real*8     :: amat_51, amat_52, amat_55, amat_56, amat_57, amat_61, amat_62, amat_65, amat_66, amat_67
real*8     :: amat_76, amat_77, amat_81, amat_85, amat_86, amat_87 
real*8     :: element_size_ij, element_size_kl, element_size_perp
real*8     :: grad_t(2), B0_R, B0_Z, factor_cs_bnd_integral, neutral_source, c_angle
logical    :: xpoint2
integer    :: n_tor_local 

type (type_node)         :: tmp_node

theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

Zbig = 1.d12

c_angle = 0.0174524d0 ! --- 1 degree angle factor for minimum heat and particle flux


!--------------------- reorder the nodes to have the same direction as full element (maybe not necesary)
if ((vertex(1) .eq. 3) .and. (vertex(2) .eq. 4)) then
  tmp_node  = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 4
  vertex(2) = 3
endif
if ((vertex(1) .eq. 4) .and. (vertex(2) .eq. 1)) then
  tmp_node  = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 1
  vertex(2) = 4
endif
if ((vertex(1) .eq. 3) .and. (vertex(2) .eq. 2)) then
  tmp_node  = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 2
  vertex(2) = 3
endif
if ((vertex(1) .eq. 2) .and. (vertex(2) .eq. 1)) then
  tmp_node  = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 1
  vertex(2) = 2
endif


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

n_tor_local = i_tor_max - i_tor_min +1
!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

  ws = wgauss(ms)

  dl   = sqrt(x_s(ms)**2 + y_s(ms)**2) 
  xjac = x_s(ms)*y_t(ms) - x_t(ms)*y_s(ms)
  BigR = x_g(ms)

  grad_t = (/ - y_s(ms),   x_s(ms) /) / xjac

!  normal_direction = (/R_mid - R_cnt, Z_mid - Z_cnt /) / norm2((/R_mid - R_cnt, Z_mid - Z_cnt /))
  normal_direction = (/x_g(ms) - R_cnt, y_g(ms) - Z_cnt /) / norm2((/x_g(ms) - R_cnt, y_g(ms) - Z_cnt /))

  normal = dot_product(grad_t,normal_direction) * grad_t      ! outward pointing normal
  normal = normal / norm2(normal)

  neutral_source = 0.d-0

  do is = 1, 10
    if     ( ((x_g(ms) - neutral_line_R_start(is))*(x_g(ms) - neutral_line_R_end(is)) .lt. 0.d0) &
       .and. ((y_g(ms) - neutral_line_Z_start(is))*(y_g(ms) - neutral_line_Z_end(is)) .lt. 0.d0) ) then
       neutral_source = neutral_source + neutral_line_source(is)
    endif
  enddo

  do mp = 1, n_plane

    ps0   = eq_g(mp,1,ms)
    ps0_s = eq_s(mp,1,ms) 
    ps0_t = eq_t(mp,1,ms)   
    ps0_x = (   y_t(ms) * ps0_s - y_s(ms) * ps0_t ) / xjac
    ps0_y = ( - x_t(ms) * ps0_s + x_s(ms) * ps0_t ) / xjac

    B0_R =   ps0_y / x_g(ms)
    B0_Z = - ps0_x / x_g(ms)

    r0    = eq_g(mp,5,ms)
    r0_s  = eq_s(mp,5,ms)
    r0_t  = eq_t(mp,5,ms)
    r0_p  = eq_p(mp,5,ms)
    r0_x = (   y_t(ms) * r0_s - y_s(ms) * r0_t ) / xjac
    r0_y = ( - x_t(ms) * r0_s + x_s(ms) * r0_t ) / xjac

    T0    = eq_g(mp,6,ms)
    T0_s  = eq_s(mp,6,ms)
    T0_t  = eq_t(mp,6,ms)
    T0_p  = eq_p(mp,6,ms)
     
    T0_x = (   y_t(ms) * T0_s - y_s(ms) * T0_t ) / xjac
    T0_y = ( - x_t(ms) * T0_s + x_s(ms) * T0_t ) / xjac

    u0_s  = eq_s(mp,2,ms)
    Vpar0 = eq_g(mp,7,ms)

    T0_corr = corr_neg_temp1(T0)
    r0_corr = corr_neg_dens(r0)

    cs0      = sqrt(gamma*T0_corr)

    Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

    BB2 = Btot**2

    bdotn = (+ ps0_y * normal(1) - ps0_x * normal(2)) / x_g(ms) / Btot

    normal_sign  = sign(1.d0,bdotn)
    normal_sign3 = sign(1.d0,ps0_s) * normal_sign

    if (vpar_smoothing) then
      factor = (0.5d0 + 0.5d0 * tanh((abs(bdotn) - 0.02d0)/0.016d0))**2 - 5.75446348E-03
    else
      factor = 1.d0
    endif

    factor_cs_bnd_integral = 0.d0
    if (mach_one_bnd_integral) factor_cs_bnd_integral = 1.d0

    do i=1,2                ! loop over nodes

      do j=1,2              ! loop over basis functions

        j2 = direction(j)
        element_size_ij = element%size(vertex(i),j2)

        do im=i_tor_min, i_tor_max

          v   =  H1(i,j,ms) * element_size_ij * HZ(im,mp)         ! test function

          rhs_ij_5 = + v * density_reflection * r0_corr * vpar0 * ps0_s * normal_sign3 * tstep  & ! right hand side equation 5
                     - v * r0_corr * cs0 * BigR * dl * c_angle * tstep                          & ! particle flux at 1 degree angle  
                     - v * r0_corr * BigR**2.d0 * u0_s * normal_sign3 * tstep                     ! reflect v_perp particle flow
           
          rhs_ij_6 = - v * (gamma_sheath -1.d0) * r0_corr * T0_corr * vpar0 * ps0_s * normal_sign3 * tstep  & ! right hand side equation 6
                     - v * (gamma_sheath -1.d0) * r0_corr * T0_corr * cs0   * BigR  * dl * c_angle * tstep  &
                     - v *                        r0_corr * T0_corr * BigR**2.d0    * u0_s  * normal_sign3 * tstep  

          rhs_ij_7 = - v * (vpar0 * Btot * normal_sign - cs0 * factor) * dl * Zbig                ! right hand side equation 7

          rhs_ij_8 = + v * neutral_reflection * r0_corr * vpar0 * ps0_s * normal_sign3 * tstep &
                     + v * neutral_reflection * r0_corr * cs0 * BigR * dl * c_angle    * tstep &  ! particle flux at 1 degree angle  
!                     - v * neutral_reflection * D_prof  * (r0_x * y_t(ms) - r0_y * x_t(ms)) * BigR * tstep &            
                     + v * neutral_source * BigR * dl * tstep                                     ! neutral source

          index_ij = n_tor_local*n_var*n_degrees*(vertex(i)-1) + n_tor_local * n_var * (j2-1) + im - i_tor_min +1  ! index in the ELM matrix

          ij5 = index_ij + 4*n_tor_local                                          ! local index in element matrix
          ij6 = index_ij + 5*n_tor_local                                          ! local index in element matrix
          ij7 = index_ij + 6*n_tor_local                                          ! local index in element matrix
          ij8 = index_ij + 7*n_tor_local                                          ! local index in element matrix

          RHS(ij5) = RHS(ij5) + rhs_ij_5 * ws                               ! add to element RHS
          RHS(ij6) = RHS(ij6) + rhs_ij_6 * ws                               ! add to element RHS
          RHS(ij7) = RHS(ij7) + rhs_ij_7 * ws * factor_cs_bnd_integral      ! add to element RHS
          RHS(ij8) = RHS(ij8) + rhs_ij_8 * ws

          do k=1,2                                                          ! loop over nodes

            do l=1,2                                                        ! loop over basis functions

              l2 = direction(l)
              element_size_kl = element%size(vertex(k),l2)

              l3 = direction_perp(j)
              element_size_perp = - element%size(vertex(k),direction_perp(1)) * 3.d0

              do in = i_tor_min, i_tor_max                                              ! loop over toroidal harmonics

                psi    = H1(k,l,ms)   * element_size_kl * HZ(in,mp)
                psi_s  = H1_s(k,l,ms) * element_size_kl * HZ(in,mp)
                psi_ss = H1_ss(k,l,ms) * element_size_kl * HZ(in,mp)
                psi_t  = H1(k,l,ms)   * element_size_kl * HZ(in,mp) * element_size_perp

                rho   = psi
                rho_s = psi_s
                rho_t = psi_t
                rho_x = (   y_t(ms) * rho_s - y_s(ms) * rho_t ) / xjac
                rho_y = ( - x_t(ms) * rho_s + x_s(ms) * rho_t ) / xjac

                T   = psi;    vpar = psi;  vpar_ss = psi_ss;  u_s = psi_s

                cs_T  = gamma * T / (2.d0 * cs0)

                amat_51 = - v * density_reflection * r0_corr  * vpar0 * psi_s * normal_sign3 * theta * tstep 

                amat_52 = + v * r0_corr * BigR**2.d0 * u_s                    * normal_sign3 * theta * tstep  

                amat_55 = - v * density_reflection * rho      * vpar0 * ps0_s * normal_sign3 * theta * tstep & 
                          + v                      * rho      * cs0   * BigR * dl * c_angle  * theta * tstep &
                          + v * rho * BigR**2.d0 * u0_s                       * normal_sign3 * theta * tstep  


                amat_56 = + v                      * r0_corr  * cs_T  * BigR * dl * c_angle  * theta * tstep
                amat_57 = - v * density_reflection * r0_corr  * vpar  * ps0_s * normal_sign3 * theta * tstep 

                amat_61 = + v * (gamma_sheath-1.d0) * r0_corr  * T0_corr * vpar0 * psi_s * normal_sign3 * theta * tstep 
                
                amat_62 = + v * r0_corr * BigR**2.d0 * u_s * normal_sign3                               * theta * tstep

                amat_65 = + v * (gamma_sheath-1.d0) * rho      * T0_corr * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                          + v * (gamma_sheath-1.d0) * rho      * T0_corr * cs0   * BigR  * dl * c_angle * theta * tstep 

                amat_66 = + v * (gamma_sheath-1.d0) * r0_corr  * T       * vpar0 * ps0_s * normal_sign3 * theta * tstep &
                          + v * (gamma_sheath-1.d0) * r0_corr  * T       * cs0   * BigR  * dl * c_angle * theta * tstep &
                          + v * (gamma_sheath-1.d0) * r0_corr  * T0_corr * cs_T  * BigR  * dl * c_angle * theta * tstep

                amat_67 = + v * (gamma_sheath-1.d0) * r0_corr  * T0_corr * vpar  * ps0_s * normal_sign3 * theta * tstep 

           
                amat_76 =   v * ( - cs_T) * factor          * dl * Zbig
                amat_77 =   v * (vpar * Btot * normal_sign) * dl * Zbig 

                amat_81 = - v * neutral_reflection * r0_corr * vpar0 * psi_s * normal_sign3      * theta * tstep 

                amat_85 = - v * neutral_reflection * rho     * vpar0 * ps0_s * normal_sign3      * theta * tstep &
                          - v * neutral_reflection * rho     * cs0         * BigR * dl * tstep * c_angle * theta * tstep 

                amat_86 = - v * neutral_reflection * r0_corr * cs_T * BigR * dl * tstep * c_angle * theta * tstep !&
 !                          + v * neutral_reflection * D_prof  * (rho_x * y_t(ms) - rho_y * x_t(ms)) * BigR * theta * tstep 

                amat_87 = - v * neutral_reflection * r0_corr * vpar  * ps0_s * normal_sign3      * theta * tstep 

                index_kl = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l2-1) + in - i_tor_min +1  ! index in the ELM matrix
                 
                kl1 = index_kl
                kl2 = index_kl + 1*n_tor_local
                kl5 = index_kl + 4*n_tor_local
                kl6 = index_kl + 5*n_tor_local
                kl7 = index_kl + 6*n_tor_local
                kl8 = index_kl + 7*n_tor_local

                ELM(ij5,kl1) =  ELM(ij5,kl1) + ws * amat_51
                ELM(ij5,kl2) =  ELM(ij5,kl2) + ws * amat_52
                ELM(ij5,kl5) =  ELM(ij5,kl5) + ws * amat_55
                ELM(ij5,kl6) =  ELM(ij5,kl6) + ws * amat_56
                ELM(ij5,kl7) =  ELM(ij5,kl7) + ws * amat_57

                ELM(ij6,kl1) =  ELM(ij6,kl1) + ws * amat_61
                ELM(ij6,kl2) =  ELM(ij6,kl2) + ws * amat_62
                ELM(ij6,kl5) =  ELM(ij6,kl5) + ws * amat_65
                ELM(ij6,kl6) =  ELM(ij6,kl6) + ws * amat_66
                ELM(ij6,kl7) =  ELM(ij6,kl7) + ws * amat_67

                ELM(ij7,kl6) =  ELM(ij7,kl6) + ws * amat_76 * factor_cs_bnd_integral
                ELM(ij7,kl7) =  ELM(ij7,kl7) + ws * amat_77 * factor_cs_bnd_integral

                ELM(ij8,kl1) =  ELM(ij8,kl1) + ws * amat_81
                ELM(ij8,kl5) =  ELM(ij8,kl5) + ws * amat_85
                ELM(ij8,kl6) =  ELM(ij8,kl6) + ws * amat_86
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
