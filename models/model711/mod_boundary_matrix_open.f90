module mod_boundary_matrix_open
  implicit none

contains

#include "corr_neg_include.f90"

subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundaries of one element
! implements the natural boundary conditions
!---------------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use diffusivities, only: get_dperp, get_zkperp

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)        ! the two nodes containing the boundary nodes
integer, intent(in)   :: i_tor_min   
integer, intent(in)   :: i_tor_max   

real*8     :: ELM(n_vertex_max*n_var*n_degrees*n_tor,n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: RHS(n_vertex_max*n_var*n_degrees*n_tor)

integer    :: vertex(2), direction(n_degrees_1d), xcase2
real*8     :: psi_axis, R_axis, Z_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
logical    :: xpoint2
real*8     :: R_g(n_gauss), R_s(n_gauss), R_t(n_gauss)
real*8     :: Z_g(n_gauss), Z_s(n_gauss), Z_t(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_t(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss), delta_t(n_plane,n_var,n_gauss)
real*8     :: Fprofile(n_gauss)

real*8     :: Qbnd(n_var), Qjac(n_var,n_var)

integer    :: i, j, j2, ms, mt, mp, k, l, l2, index_ij, index_kl, ij, kl
integer    :: in, im, ivar, kvar
integer    :: j3, direction_perp(n_degrees_1d)
real*8     :: ws, xjac,  R, phi, DL, Zbig
real*8     :: R_mid, Z_mid, R_cnt, Z_cnt
real*8     :: theta, zeta, psi_norm, ZK_prof, integrand

real*8     :: c_s, cs_Ti, cs_Te
real*8     :: AR0, AR0_p, AR0_s, AR0_t, AR0_R, AR0_Z
real*8     :: AZ0, AZ0_p, AZ0_s, AZ0_t, AZ0_R, AZ0_Z     
real*8     :: A30, A30_p, A30_s, A30_t, A30_R, A30_Z
real*8     :: uR0, uR0_s, uR0_t, uR0_R, uR0_Z
real*8     :: uZ0, uZ0_s, uZ0_t, uZ0_R, uZ0_Z
real*8     :: UP0, Up
real*8     :: rho0,rho0_p,rho0_s,rho0_t,rho0_R,rho0_Z
real*8     :: Ti0, Ti0_p, Ti0_s, Ti0_t, Ti0_R, Ti0_Z, Ti0_corr
real*8     :: Te0, Te0_p, Te0_s, Te0_t, Te0_R, Te0_Z, Te0_corr
real*8     :: pi0,               pi0_p, pi0_R, pi0_Z
real*8     :: AR, AR_p, AR_s, AR_t, AR_R, AR_Z
real*8     :: AZ, AZ_p, AZ_s, AZ_t, AZ_R, AZ_Z     
real*8     :: A3, A3_p, A3_s, A3_t, A3_R, A3_Z
real*8     :: uR, uR_s, uR_t, uR_R, uR_Z
real*8     :: uZ, uZ_s, uZ_t, uZ_R, uZ_Z
real*8     :: rho,rho_p,rho_s,rho_t,rho_R,rho_Z
real*8     :: Ti, Ti_p, Ti_s, Ti_t, Ti_R, Ti_Z
real*8     :: Te, Te_p, Te_s, Te_t, Te_R, Te_Z
real*8     :: bf, bf_s, bf_t, bf_p, bf_R, bf_Z

real*8     :: BB2, BB2_AR, BB2_AZ, BB2_A3
real*8     :: BR0, BR0_AR, BR0_AZ, BR0_A3
real*8     :: BZ0, BZ0_AR, BZ0_AZ, BZ0_A3
real*8     :: Bp0, Bp0_AR, Bp0_AZ, Bp0_A3

real*8     :: B_dot_n, B_dot_n_AR, B_dot_n_AZ, B_dot_n_A3, cs_direction
real*8     :: ZKpar_Ti, ZKpar_Te

real*8     :: rhoVdiaR0, rhoVdiaR0_AR, rhoVdiaR0_AZ, rhoVdiaR0_A3, rhoVdiaR0_rho, rhoVdiaR0_Ti
real*8     :: rhoVdiaZ0, rhoVdiaZ0_AR, rhoVdiaZ0_AZ, rhoVdiaZ0_A3, rhoVdiaZ0_rho, rhoVdiaZ0_Ti
real*8     :: rhoVdia_dot_n, rhoVdia_dot_n_AR, rhoVdia_dot_n_AZ, rhoVdia_dot_n_A3, rhoVdia_dot_n_rho, rhoVdia_dot_n_Ti

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: element_size_ij, element_size_kl, element_size_perp
real*8     :: normal(2), normal_direction(2)
real*8     :: grad_s(2), grad_t(2)
real*8     :: Mach1
integer    :: n_tor_local

type (type_node)         :: tmp_node

logical    :: parallel_projection

! --- Time integration parameters
theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

! --- Needs adaptation for t-derivatives
if (direction(2) == 3) return

! --- Flag to switch Mach-1 between boundary_conditions and boundary_matrix_open
Mach1 = 0.d0
if (Mach1_openBC) Mach1 = 1.d0
parallel_projection = .false. !.true. ! note this is not exactly the same as the projection of the momentum equation, so we keep the option here...

! --- Penalisation cofficient to impose BCs
zbig = 1.d11

! --- Initialise variables before integration
R_g     = 0.d0; R_s     = 0.d0; R_t     = 0.d0; 
Z_g     = 0.d0; Z_s     = 0.d0; Z_t     = 0.d0;
eq_g    = 0.d0; eq_s    = 0.d0; eq_t    = 0.d0; eq_p = 0.d0;
delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0;
Fprofile = 0.d0

! --- Strategic points on elements to define normal vector properly
R_mid = sum(nodes(1:2)%x(1,1,1)) / 2.d0     ! mid point on boundary (approx.)
Z_mid = sum(nodes(1:2)%x(1,1,2)) / 2.d0
R_cnt = sum(nodes(1:4)%x(1,1,1)) / 4.d0     ! center point within element (approx.)
Z_cnt = sum(nodes(1:4)%x(1,1,2)) / 4.d0

normal_direction = (/R_mid - R_cnt, Z_mid - Z_cnt /) / norm2((/R_mid - R_cnt, Z_mid - Z_cnt /))
direction_perp(1) = 6 / direction(2)     ! =3 if direction(2)=2, =3 if direction(2)=3
direction_perp(2) = 4
if (n_order .ge. 5) then
  direction_perp(1) = 6 / direction(2)     ! =3 if direction(2)=2, =3 if direction(2)=3
  direction_perp(2) = 4
  if (direction(2) .eq. 2) direction_perp(3) = 7
  if (direction(2) .eq. 3) direction_perp(3) = 8
endif

! --- Loop over nodes
do i=1,2
  ! --- Loop over basis functions
  do j=1,n_degrees_1d

    j2 = direction(j)
    j3 = direction_perp(j)

    element_size_ij = element%size(vertex(i),j2)
    !element_size_perp = - element%size(vertex(i),direction_perp(1)) * 3.d0
    if(vertex(1) == 1)then ! NEEDS TO BE CONFIRMED FOR WALL-GRIDS AND t-derivaties !
      element_size_perp = + element%size(vertex(i),j3) * 3.d0
    elseif(vertex(1)==3)then
      element_size_perp = - element%size(vertex(i),j3) * 3.d0
    endif

    ! --- Gaussian integration
    do ms=1, n_gauss

      ! --- Pre-define R and Z
      R_g(ms)  = R_g(ms)  + nodes(i)%x(1,j2,1) * element_size_ij * H1(i,j,ms)
      R_s(ms)  = R_s(ms)  + nodes(i)%x(1,j2,1) * element_size_ij * H1_s(i,j,ms)
      R_t(ms)  = R_t(ms)  + nodes(i)%x(1,j3,1) * element_size_ij * H1(i,j,ms)   * element_size_perp

      Z_g(ms)  = Z_g(ms)  + nodes(i)%x(1,j2,2) * element_size_ij * H1(i,j,ms)
      Z_s(ms)  = Z_s(ms)  + nodes(i)%x(1,j2,2) * element_size_ij * H1_s(i,j,ms)
      Z_t(ms)  = Z_t(ms)  + nodes(i)%x(1,j3,2) * element_size_ij * H1(i,j,ms)   * element_size_perp

      ! --- Pre-define F-profile from initial equilibrium
      Fprofile(ms)   = Fprofile(ms)   + nodes(i)%Fprof_eq(j2)    * element_size_ij * H1(i,j,ms)

      ! --- Toroidal integration for each variable
      do mp=1,n_plane
        do k=1,n_var
          do in=1,n_tor

            ! --- Save variables
            eq_g(mp,k,ms)  = eq_g(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ  (in,mp)
            eq_s(mp,k,ms)  = eq_s(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ  (in,mp)
            eq_t(mp,k,ms)  = eq_t(mp,k,ms)  + nodes(i)%values(in,j3,k) * element_size_ij * H1(i,j,ms)   * HZ  (in,mp) * element_size_perp
            eq_p(mp,k,ms)  = eq_p(mp,k,ms)  + nodes(i)%values(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ_p(in,mp)

            ! --- Save deltas
            delta_g(mp,k,ms) = delta_g(mp,k,ms) + nodes(i)%deltas(in,j2,k) * element_size_ij * H1(i,j,ms)   * HZ(in,mp)
            delta_s(mp,k,ms) = delta_s(mp,k,ms) + nodes(i)%deltas(in,j2,k) * element_size_ij * H1_s(i,j,ms) * HZ(in,mp)

          enddo
        enddo
      enddo

    enddo
  enddo
enddo

n_tor_local = i_tor_max - i_tor_min + 1
! --- Gaussian integration
do ms=1, n_gauss

  ! --- Gaussian weight
  ws = wgauss(ms)

  ! --- Jacobian
  xjac = R_s(ms) * Z_t(ms) - R_t(ms) * Z_s(ms)

  ! --- Normal at target
  grad_s = (/   Z_t(ms), - R_t(ms) /) / xjac
  grad_t = (/ - Z_s(ms),   R_s(ms) /) / xjac
  normal = dot_product(grad_t,normal_direction) * grad_t
  normal = normal / norm2(normal)

  ! --- Curve integrand
  DL   = sqrt(R_s(ms)**2 + Z_s(ms)**2)
  R = R_g(ms)
  integrand = ws * R * DL ! IMPORTANT: is DL ignored in model303 ?

  ! --- Intitilise RHS and LHS
  Qbnd = 0.d0
  Qjac = 0.d0

  ! --- Toroidal integration
  do mp = 1, n_plane

    ! --- Density
    ! --- Density
    rho0    = eq_g(mp,var_rho,ms)
    rho0_p  = eq_p(mp,var_rho,ms)
    rho0_s  = eq_s(mp,var_rho,ms)
    rho0_t  = eq_t(mp,var_rho,ms)
    rho0_R  = (   Z_t(ms) * rho0_s  - Z_s(ms) * rho0_t ) / xjac
    rho0_Z  = ( - R_t(ms) * rho0_s  + R_s(ms) * rho0_t ) / xjac
    
    ! --- Velocity
    uR0   = eq_g(mp,var_uR,ms)
    uZ0   = eq_g(mp,var_uZ,ms)
    up0   = eq_g(mp,var_up,ms)

    ! --- Ti
    Ti0    = eq_g(mp,var_Ti,ms)
    Ti0_p  = eq_p(mp,var_Ti,ms)
    Ti0_s  = eq_s(mp,var_Ti,ms)
    Ti0_t  = eq_t(mp,var_Ti,ms)
    Ti0_R  = (   Z_t(ms) * Ti0_s  - Z_s(ms) * Ti0_t ) / xjac
    Ti0_Z  = ( - R_t(ms) * Ti0_s  + R_s(ms) * Ti0_t ) / xjac

    ! --- Te
    Te0    = eq_g(mp,var_Te,ms)
    Te0_p  = eq_p(mp,var_Te,ms)
    Te0_s  = eq_s(mp,var_Te,ms)
    Te0_t  = eq_t(mp,var_Te,ms)
    Te0_R  = (   Z_t(ms) * Te0_s  - Z_s(ms) * Te0_t ) / xjac
    Te0_Z  = ( - R_t(ms) * Te0_s  + R_s(ms) * Te0_t ) / xjac

    ! --- P
    pi0    = rho0 * Ti0
    pi0_R  = rho0_R * Ti0 + rho0 * Ti0_R
    pi0_Z  = rho0_Z * Ti0 + rho0 * Ti0_Z
    pi0_p  = rho0_p * Ti0 + rho0 * Ti0_p

    ! --- AR
    AR0   = eq_g(mp,var_AR,ms)
    AR0_p = eq_p(mp,var_AR,ms)
    AR0_s = eq_s(mp,var_AR,ms)
    AR0_t = eq_t(mp,var_AR,ms)
    AR0_R = (   Z_t(ms) * AR0_s  - Z_s(ms) * AR0_t ) / xjac
    AR0_Z = ( - R_t(ms) * AR0_s  + R_s(ms) * AR0_t ) / xjac

    ! --- AZ
    AZ0   = eq_g(mp,var_AZ,ms)
    AZ0_p = eq_p(mp,var_AZ,ms)
    AZ0_s = eq_s(mp,var_AZ,ms)
    AZ0_t = eq_t(mp,var_AZ,ms)
    AZ0_R = (   Z_t(ms) * AZ0_s  - Z_s(ms) * AZ0_t ) / xjac
    AZ0_Z = ( - R_t(ms) * AZ0_s  + R_s(ms) * AZ0_t ) / xjac

    ! --- A3
    A30   = eq_g(mp,var_A3,ms)
    A30_p = eq_p(mp,var_A3,ms)
    A30_s = eq_s(mp,var_A3,ms)
    A30_t = eq_t(mp,var_A3,ms)
    A30_R = (   Z_t(ms) * A30_s  - Z_s(ms) * A30_t ) / xjac
    A30_Z = ( - R_t(ms) * A30_s  + R_s(ms) * A30_t ) / xjac

    ! --- Magnetic field
    BR0 = ( A30_Z - AZ0_p )/ R
    BZ0 = ( AR0_p - A30_R )/ R
    Bp0 = ( AZ0_R - AR0_Z )       +   Fprofile(ms) / R

    BB2 = BR0*BR0 + BZ0*BZ0 + Bp0*Bp0

    ! --- Diamagnetic velocity (times rho)
    rhoVdiaR0 = tauIC*F0 / (R * BB2) * (  BZ0*pi0_p - R*Bp0*pi0_Z)
    rhoVdiaZ0 = tauIC*F0 / (R * BB2) * (R*BP0*pi0_R -   BR0*pi0_p)
    rhoVdia_dot_n = rhoVdiaR0 * normal(1) + rhoVdiaZ0 * normal(2)

    ! --- Magnetic field direction at target
    B_dot_n = BR0 * normal(1) + BZ0 * normal(2)
    cs_direction = B_dot_n / abs(B_dot_n)

    Ti0_corr = max(Ti0,1.d-12) ! CAREFUL! FULL-MHD DOESN'T LIKE THE CORR FUNCTIONS AT ALL
    Te0_corr = max(Te0,1.d-12) ! CAREFUL! FULL-MHD DOESN'T LIKE THE CORR FUNCTIONS AT ALL
    c_s = sqrt(gamma * (Ti0_corr+Te0_corr))

    ! --- Loop over nodes
    do i=1,2
      ! --- Loop over basis functions
      do j=1,n_degrees_1d

        j2 = direction(j)
        element_size_ij = element%size(vertex(i),j2)

        ! Loop over toroidal modes
        do im=i_tor_min, i_tor_max

          ! --- Test function
          v   =  H1(i,j,ms) * element_size_ij * HZ(im,mp)

          ! --- VR-equation
          Qbnd(var_uR)   = Mach1 * zbig * v * ( UR0  - c_s * BR0 * cs_direction / sqrt(BB2) )

          ! --- VZ-equation
          Qbnd(var_uZ)   = Mach1 * zbig * v * ( UZ0  - c_s * BZ0 * cs_direction / sqrt(BB2) )

          ! --- Vp-equation
          if (parallel_projection) then
            Qbnd(var_up) = Mach1 * zbig * v * ( BR0 * UR0 + BZ0 * UZ0 + Bp0 * Up0 - c_s * cs_direction * sqrt(BB2) )
          else
            Qbnd(var_up) = Mach1 * zbig * v * ( Up0  - c_s * Bp0 * cs_direction / sqrt(BB2) )
          endif

          ! --- Diamagnetic BC's
          Qbnd(var_rho) = - v * rhoVdia_dot_n

          ! --- Sheath BC's
          Qbnd(var_Ti) = - v * (gamma_sheath_i - 1.d0) * rho0 * Ti0 * cs_direction * c_s * B_dot_n / sqrt(BB2)
          Qbnd(var_Te) = - v * (gamma_sheath_e - 1.d0) * rho0 * Te0 * cs_direction * c_s * B_dot_n / sqrt(BB2)

          ! --- Fill in RHS
          index_ij = n_tor_local*n_var*n_degrees*(vertex(i)-1) + n_tor_local * n_var * (j2-1) + im - i_tor_min +1  ! index in the ELM matrix
          do ivar= 1,n_var
            ij = index_ij + (ivar-1)*n_tor_local
            RHS(ij) =  RHS(ij) + Qbnd(ivar) * integrand * tstep
          enddo

          ! --- loop over nodes
          do k=1,2

            ! --- loop over basis functions
            do l=1,n_degrees_1d

              l2 = direction(l)

              element_size_kl   = element%size(vertex(k),l2)
              !element_size_perp = - element%size(vertex(k),direction_perp(1)) * 3.d0
              if(vertex(1) == 1)then ! NEEDS TO BE CONFIRMED FOR WALL-GRIDS AND t-derivaties !
                element_size_perp = + element%size(vertex(k),direction_perp(l)) * 3.d0
              elseif(vertex(1)==3)then
                element_size_perp = - element%size(vertex(k),direction_perp(l)) * 3.d0
              endif

              ! Loop over toroidal modes
              do in=i_tor_min, i_tor_max

                ! --- Basis functions
                bf   = H1(k,l,ms)   * element_size_kl * HZ(in,mp)
                bf_s = H1_s(k,l,ms) * element_size_kl * HZ(in,mp)   
                bf_t = H1(k,l,ms)   * element_size_kl * HZ(in,mp) * element_size_perp
                bf_p = H1(k,l,ms)   * element_size_kl * HZ_p(in,mp)
                bf_R = (   Z_t(ms) * bf_s - Z_s(ms) * bf_t ) / xjac
                bf_Z = ( - R_t(ms) * bf_s + R_s(ms) * bf_t ) / xjac

                ! --- Copies of basis functions
                uR    = bf    ;  uZ    = bf    ;  up    = bf
                AR    = bf    ;  AZ    = bf    ;  A3    = bf    ;  rho   = bf    ;  Ti   = bf    ;  Te   = bf
                AR_R  = bf_R  ;  AZ_R  = bf_R  ;  A3_R  = bf_R  ;  rho_R = bf_R  ;  Ti_R = bf_R  ;  Te_R = bf_R
                AR_Z  = bf_Z  ;  AZ_Z  = bf_Z  ;  A3_Z  = bf_Z  ;  rho_Z = bf_Z  ;  Ti_Z = bf_Z  ;  Te_Z = bf_Z
                AR_p  = bf_p  ;  AZ_p  = bf_p  ;  A3_p  = bf_p  ;  rho_p = bf_p  ;  Ti_p = bf_p  ;  Te_p = bf_p
                AR_s  = bf_s  ;  AZ_s  = bf_s  ;  A3_s  = bf_s  ;  rho_s = bf_s  ;  Ti_s = bf_s  ;  Te_s = bf_s
                AR_t  = bf_t  ;  AZ_t  = bf_t  ;  A3_t  = bf_t  ;  rho_t = bf_t  ;  Ti_t = bf_t  ;  Te_t = bf_t

                ! --- Magnetic field derivatives
                BR0_AR =   0.d0     ; BR0_AZ = - AZ_p / R ; BR0_A3 =   A3_Z / R
                BZ0_AR =   AR_p / R ; BZ0_AZ =   0.d0     ; BZ0_A3 = - A3_R / R
                Bp0_AR = - AR_Z     ; Bp0_AZ =   AZ_R     ; Bp0_A3 =   0.d0

                B_dot_n_AR = BR0_AR * normal(1) + BZ0_AR * normal(2)
                B_dot_n_AZ = BR0_AZ * normal(1) + BZ0_AZ * normal(2)
                B_dot_n_A3 = BR0_A3 * normal(1) + BZ0_A3 * normal(2)

                BB2_AR = 2.d0*(BR0_AR * BR0 + BZ0_AR * BZ0 + Bp0_AR * Bp0 )
                BB2_AZ = 2.d0*(BR0_AZ * BR0 + BZ0_AZ * BZ0 + Bp0_AZ * Bp0 )
                BB2_A3 = 2.d0*(BR0_A3 * BR0 + BZ0_A3 * BZ0 + Bp0_A3 * Bp0 )

                ! --- Diamagnetic BC's variables
                rhoVdiaR0_AR  = + tauIC*F0 / (R * BB2   ) * (  BZ0_AR*pi0_p - R*Bp0_AR*pi0_Z) &
                                - tauIC*F0 / (R * BB2**2) * (  BZ0   *pi0_p - R*Bp0   *pi0_Z) * BB2_AR
                rhoVdiaR0_AZ  = + tauIC*F0 / (R * BB2   ) * (  BZ0_AZ*pi0_p - R*Bp0_AZ*pi0_Z) &
                                - tauIC*F0 / (R * BB2**2) * (  BZ0   *pi0_p - R*Bp0   *pi0_Z) * BB2_AZ
                rhoVdiaR0_A3  = + tauIC*F0 / (R * BB2   ) * (  BZ0_A3*pi0_p - R*Bp0_A3*pi0_Z) &
                                - tauIC*F0 / (R * BB2**2) * (  BZ0   *pi0_p - R*Bp0   *pi0_Z) * BB2_A3
                rhoVdiaR0_rho = + tauIC*F0 / (R * BB2   ) * (  BZ0*(rho*Ti0_p+rho_p*Ti0) - R*Bp0*(rho*Ti0_Z+rho_Z*Ti0))
                rhoVdiaR0_Ti  = + tauIC*F0 / (R * BB2   ) * (  BZ0*(rho0*Ti_p+rho0_p*Ti) - R*Bp0*(rho0*Ti_Z+rho0_Z*Ti))

                rhoVdiaZ0_AR  = + tauIC*F0 / (R * BB2   ) * (R*BP0_AR*pi0_R -   BR0_AR*pi0_p) &
                                - tauIC*F0 / (R * BB2**2) * (R*BP0   *pi0_R -   BR0   *pi0_p) * BB2_AR
                rhoVdiaZ0_AZ  = + tauIC*F0 / (R * BB2   ) * (R*BP0_AZ*pi0_R -   BR0_AZ*pi0_p) &
                                - tauIC*F0 / (R * BB2**2) * (R*BP0   *pi0_R -   BR0   *pi0_p) * BB2_AZ
                rhoVdiaZ0_A3  = + tauIC*F0 / (R * BB2   ) * (R*BP0_A3*pi0_R -   BR0_A3*pi0_p) &
                                - tauIC*F0 / (R * BB2**2) * (R*BP0   *pi0_R -   BR0   *pi0_p) * BB2_A3
                rhoVdiaZ0_rho = + tauIC*F0 / (R * BB2   ) * (R*BP0*(rho*Ti0_R+rho_R*Ti0) -   BR0*(rho*Ti0_p+rho_p*Ti0))
                rhoVdiaZ0_Ti  = + tauIC*F0 / (R * BB2   ) * (R*BP0*(rho0*Ti_R+rho0_R*Ti) -   BR0*(rho0*Ti_p+rho0_p*Ti))

                rhoVdia_dot_n_AR  = rhoVdiaR0_AR  * normal(1) + rhoVdiaZ0_AR  * normal(2)
                rhoVdia_dot_n_AZ  = rhoVdiaR0_AZ  * normal(1) + rhoVdiaZ0_AZ  * normal(2)
                rhoVdia_dot_n_A3  = rhoVdiaR0_A3  * normal(1) + rhoVdiaZ0_A3  * normal(2)
                rhoVdia_dot_n_rho = rhoVdiaR0_rho * normal(1) + rhoVdiaZ0_rho * normal(2)
                rhoVdia_dot_n_Ti  = rhoVdiaR0_Ti  * normal(1) + rhoVdiaZ0_Ti  * normal(2)

                ! --- Mach-1 BC's
                cs_Ti = gamma * Ti / (2.d0 * c_s)
                cs_Te = gamma * Te / (2.d0 * c_s)

                ! --- VR-linearised equation
                Qjac(var_uR,var_uR) = - Mach1 * zbig * v * UR
                Qjac(var_uR,var_AR) = - Mach1 * zbig * v * c_s * cs_direction * ( - BR0_AR / sqrt(BB2) + 0.5 * BR0 * BB2_AR / BB2**1.5 )
                Qjac(var_uR,var_AZ) = - Mach1 * zbig * v * c_s * cs_direction * ( - BR0_AZ / sqrt(BB2) + 0.5 * BR0 * BB2_AZ / BB2**1.5 )
                Qjac(var_uR,var_A3) = - Mach1 * zbig * v * c_s * cs_direction * ( - BR0_A3 / sqrt(BB2) + 0.5 * BR0 * BB2_A3 / BB2**1.5 )
                Qjac(var_uR,var_Ti) = - Mach1 * zbig * v * ( - cs_Ti * BR0 * cs_direction / sqrt(BB2) )
                Qjac(var_uR,var_Te) = - Mach1 * zbig * v * ( - cs_Te * BR0 * cs_direction / sqrt(BB2) )

                ! --- VZ-linearised equation
                Qjac(var_uZ,var_uZ) = - Mach1 * zbig * v * UZ
                Qjac(var_uZ,var_AR) = - Mach1 * zbig * v * c_s * cs_direction * ( - BZ0_AR / sqrt(BB2) + 0.5 * BZ0 * BB2_AR / BB2**1.5 )
                Qjac(var_uZ,var_AZ) = - Mach1 * zbig * v * c_s * cs_direction * ( - BZ0_AZ / sqrt(BB2) + 0.5 * BZ0 * BB2_AZ / BB2**1.5 )
                Qjac(var_uZ,var_A3) = - Mach1 * zbig * v * c_s * cs_direction * ( - BZ0_A3 / sqrt(BB2) + 0.5 * BZ0 * BB2_A3 / BB2**1.5 )
                Qjac(var_uZ,var_Ti) = - Mach1 * zbig * v * ( - cs_Ti * BZ0 * cs_direction / sqrt(BB2) )
                Qjac(var_uZ,var_Te) = - Mach1 * zbig * v * ( - cs_Te * BZ0 * cs_direction / sqrt(BB2) )

                ! --- Vp-linearised equation
                if (parallel_projection) then
                  Qjac(var_up,var_uR) = - Mach1 * zbig * v * BR0 * UR
                  Qjac(var_up,var_uZ) = - Mach1 * zbig * v * BZ0 * UZ
                  Qjac(var_up,var_up) = - Mach1 * zbig * v * Bp0 * Up
                  Qjac(var_up,var_AR) = - Mach1 * zbig * v * ( - c_s  * cs_direction * 0.5 * BB2_AR / sqrt(BB2) )
                  Qjac(var_up,var_AZ) = - Mach1 * zbig * v * ( - c_s  * cs_direction * 0.5 * BB2_AZ / sqrt(BB2) )
                  Qjac(var_up,var_A3) = - Mach1 * zbig * v * ( - c_s  * cs_direction * 0.5 * BB2_A3 / sqrt(BB2) )
                  Qjac(var_up,var_Ti) = - Mach1 * zbig * v * ( - cs_Ti * cs_direction * sqrt(BB2) )
                  Qjac(var_up,var_Te) = - Mach1 * zbig * v * ( - cs_Te * cs_direction * sqrt(BB2) )
                else
                  Qjac(var_up,var_up) = - Mach1 * zbig * v * Up
                  Qjac(var_up,var_AR) = - Mach1 * zbig * v * c_s * cs_direction * ( - Bp0_AR / sqrt(BB2) + 0.5 * Bp0 * BB2_AR / BB2**1.5 )
                  Qjac(var_up,var_AZ) = - Mach1 * zbig * v * c_s * cs_direction * ( - Bp0_AZ / sqrt(BB2) + 0.5 * Bp0 * BB2_AZ / BB2**1.5 )
                  Qjac(var_up,var_A3) = - Mach1 * zbig * v * c_s * cs_direction * ( - Bp0_A3 / sqrt(BB2) + 0.5 * Bp0 * BB2_A3 / BB2**1.5 )
                  Qjac(var_up,var_Ti) = - Mach1 * zbig * v * ( - cs_Ti * Bp0 * cs_direction / sqrt(BB2) )
                  Qjac(var_up,var_Te) = - Mach1 * zbig * v * ( - cs_Te * Bp0 * cs_direction / sqrt(BB2) )
                endif

                ! --- Diamagnetic BC's
                Qjac(var_rho, var_AR ) = + v * rhoVdia_dot_n_AR 
                Qjac(var_rho, var_AZ ) = + v * rhoVdia_dot_n_AZ 
                Qjac(var_rho, var_A3 ) = + v * rhoVdia_dot_n_A3 
                Qjac(var_rho, var_rho) = + v * rhoVdia_dot_n_rho
                Qjac(var_rho, var_Ti ) = + v * rhoVdia_dot_n_Ti 

                ! --- Sheath-linearised energy equation
                Qjac(var_Ti, var_AR )  = + v * (gamma_sheath - 1.d0) * rho0 * Ti0 * cs_direction * c_s   * B_dot_n_AR / sqrt(BB2) &
                                         - v * (gamma_sheath - 1.d0) * rho0 * Ti0 * cs_direction * c_s   * B_dot_n    * 0.5 * BB2_AR / BB2**1.5
                Qjac(var_Ti, var_AZ )  = + v * (gamma_sheath - 1.d0) * rho0 * Ti0 * cs_direction * c_s   * B_dot_n_AZ / sqrt(BB2) &
                                         - v * (gamma_sheath - 1.d0) * rho0 * Ti0 * cs_direction * c_s   * B_dot_n    * 0.5 * BB2_AZ / BB2**1.5
                Qjac(var_Ti, var_A3 )  = + v * (gamma_sheath - 1.d0) * rho0 * Ti0 * cs_direction * c_s   * B_dot_n_A3 / sqrt(BB2) &
                                         - v * (gamma_sheath - 1.d0) * rho0 * Ti0 * cs_direction * c_s   * B_dot_n    * 0.5 * BB2_A3 / BB2**1.5
                Qjac(var_Ti, var_rho)  = + v * (gamma_sheath - 1.d0) * rho  * Ti0 * cs_direction * c_s   * B_dot_n    / sqrt(BB2)
                Qjac(var_Ti, var_Ti)   = + v * (gamma_sheath - 1.d0) * rho0 * Ti  * cs_direction * c_s   * B_dot_n    / sqrt(BB2) &
                                         + v * (gamma_sheath - 1.d0) * rho0 * Ti0 * cs_direction * cs_Ti * B_dot_n    / sqrt(BB2)
                Qjac(var_Ti, var_Te)   = + v * (gamma_sheath - 1.d0) * rho0 * Ti0 * cs_direction * cs_Te * B_dot_n    / sqrt(BB2)

                ! --- Sheath-linearised energy equation
                Qjac(var_Te, var_AR )  = + v * (gamma_sheath - 1.d0) * rho0 * Te0 * cs_direction * c_s   * B_dot_n_AR / sqrt(BB2) &
                                         - v * (gamma_sheath - 1.d0) * rho0 * Te0 * cs_direction * c_s   * B_dot_n    * 0.5 * BB2_AR / BB2**1.5
                Qjac(var_Te, var_AZ )  = + v * (gamma_sheath - 1.d0) * rho0 * Te0 * cs_direction * c_s   * B_dot_n_AZ / sqrt(BB2) &
                                         - v * (gamma_sheath - 1.d0) * rho0 * Te0 * cs_direction * c_s   * B_dot_n    * 0.5 * BB2_AZ / BB2**1.5
                Qjac(var_Te, var_A3 )  = + v * (gamma_sheath - 1.d0) * rho0 * Te0 * cs_direction * c_s   * B_dot_n_A3 / sqrt(BB2) &
                                         - v * (gamma_sheath - 1.d0) * rho0 * Te0 * cs_direction * c_s   * B_dot_n    * 0.5 * BB2_A3 / BB2**1.5
                Qjac(var_Te, var_rho)  = + v * (gamma_sheath - 1.d0) * rho  * Te0 * cs_direction * c_s   * B_dot_n    / sqrt(BB2)
                Qjac(var_Te, var_Te)   = + v * (gamma_sheath - 1.d0) * rho0 * Te  * cs_direction * c_s   * B_dot_n    / sqrt(BB2) &
                                         + v * (gamma_sheath - 1.d0) * rho0 * Te0 * cs_direction * cs_Te * B_dot_n    / sqrt(BB2)
                Qjac(var_Te, var_Ti)   = + v * (gamma_sheath - 1.d0) * rho0 * Te0 * cs_direction * cs_Ti * B_dot_n    / sqrt(BB2)

                ! --- Fill-in Matrix
                index_kl = n_tor_local*n_var*n_degrees*(vertex(k)-1) + n_tor_local * n_var * (l2-1) + in - i_tor_min +1! index in the ELM matrix 
                do ivar= 1,n_var
                  do kvar= 1,n_var
                    ij = index_ij + (ivar-1)*n_tor_local
                    kl = index_kl + (kvar-1)*n_tor_local
                    ELM(ij,kl) =  ELM(ij,kl) + Qjac(ivar,kvar) * integrand * theta * tstep
                  enddo
                enddo

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


