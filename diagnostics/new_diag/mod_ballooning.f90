!> Infinite-n ideal MHD ballooning stability of the axisymmetric (n=0) part of a JOREK state.
!!
!! Ballooning equation (Connor, Hastie, Taylor, Proc. R. Soc. A 365, 1 (1979)) along a field line
!! extended in the straight field line angle theta (= theta*), with k_perp = grad(alpha):
!!
!!   d/dtheta [ g dX/dtheta ] + c X = lambda w X,      lambda = -omega^2 = gamma^2
!!
!!   g = |grad alpha|^2 / (J B^2)
!!   c = 2 J (grad alpha x B . kappa) (grad alpha x B . grad psi) P0'(psi) / B^4
!!   w = J rho |grad alpha|^2 / B^2
!!
!! with J = q R^2 / F the straight field line Jacobian, P0 = mu0 p and rho the normalized density
!! (JOREK units, lambda in units of 1/t0^2). The drive is written as a product of two terms so
!! that it does not depend on the signs of psi, F or of the poloidal orientation.
!!
!! Field line label alpha = phi - nu(psi,l), with l the poloidal arc length along the flux surface
!! and orthogonal coordinates (psi, l), i.e. grad(psi).grad(l) = 0 (as in Forest, Classic Problems
!! in MHD, lecture 26, Eqs. 26.6-26.11, where nu = mu0 I J / R^2). Then
!!
!!   grad(alpha) = grad(phi) - nu_l t - nu_psi grad(psi),  nu_l = B.grad(phi) / B.t
!!   |grad alpha|^2 = B^2/|grad psi|^2 + nu_psi^2 |grad psi|^2
!!
!! where t is the unit tangent in the direction of increasing theta* and nu_psi = d nu / d psi at
!! fixed l (the integrated local shear). Along the normal n = grad(psi)/|grad psi| the arc length
!! element of the neighbouring surface grows with the contour curvature kappa_c = t.H.t/|grad psi|
!! (H = Hessian of psi), which gives
!!
!!   nu_psi(theta) = int_theta0^theta q/|grad psi| [ n.grad(F)/F - n_R/R - n.H.n/|grad psi|
!!                                                     + t.H.t/|grad psi| ] dtheta*
!!
!! using dl/dtheta* = q/nu_l. Over one poloidal turn nu_psi increases by exactly 2 pi dq/dpsi; this
!! is used as a consistency check and gives the magnetic shear without radial finite differences.
!! The integration constant (nu_psi = 0 at theta0) is the ballooning angle theta0.
!!
!! Pressure-gradient scans (alpha_crit) use a local equilibrium variation in the spirit of
!! Greene & Chance, Nucl. Fusion 21, 453 (1981): the surface shape, |grad psi| and F on the surface
!! are kept, P0' -> P0' + dP and FF' -> FF' + dFF, which by Grad-Shafranov change
!! n.H.n by -(R^2 dP + dFF) and n.grad(F)/F by dFF |grad psi|/F^2. dFF is chosen to keep q' (or to
!! obtain a prescribed shear). This reproduces the alpha sin(theta) term of the s-alpha model and
!! thereby second stability. The actual state (dP = dFF = 0) uses the JOREK fields unchanged.
!!
!! Discretisation as in lecture 26 (Eq. 26.86): conservative second order finite differences on
!! theta0 + [-2 pi nturns, 2 pi nturns] with X = 0 at both ends. The generalized problem with
!! diagonal w is symmetrised, and eigenvalues are located with Sturm sequence counts (number of
!! unstable modes = number of eigenvalues > 0), bisection for lambda_max, inverse iteration for X.
module mod_ballooning

  use constants
  use mod_parameters
  use data_structure
  use equil_info
  use mod_new_diag

  implicit none

  private
  public :: t_ball_surface, t_ball_case, ball_surfaces_from_state, ball_evaluate, ball_alpha_crit,  &
    ball_eigenfunction, ball_salpha_lambda, ball_alpha_to_P0p, ball_cleanup_surface

  character(len=16), parameter :: THIS_MOD_NAME = 'mod_ballooning'

  !> Local quantities on one flux surface, at np points equidistant in theta* (index 0..np-1)
  type t_ball_surface
    logical :: valid = .false.
    integer :: np = 0
    real*8  :: psin = 0.d0, q = 0.d0       !< q is signed: sign of dphi/dtheta* along B
    real*8  :: volume = 0.d0, dV_dpsi = 0.d0, r_minor = 0.d0
    real*8  :: P0p = 0.d0                  !< Surface average of grad(P0).grad(psi)/|grad psi|^2
    real*8  :: P0p_var = 0.d0              !< RMS poloidal variation of that quantity (diagnostic)
    real*8  :: FFp = 0.d0                  !< Surface average of FF' = -Delta*psi - R^2 P0' (GS)
    real*8  :: FFp_var = 0.d0              !< RMS poloidal variation of FF', relative to |FF'| +
                                           !< <R^2>|P0'| (zero in a Grad-Shafranov equilibrium)
    real*8  :: alpha0 = 0.d0               !< alpha_miller of the actual state
    real*8  :: dq_dpsi = 0.d0              !< From the closure of nu_psi (= I0(np)/(2 pi))
    real*8  :: dq_dpsi_fd = 0.d0           !< From finite differences of q between surfaces
    real*8, allocatable :: theta(:), R(:), Z(:), Jac(:), BB2(:), gpsi2(:), rho(:), D0(:), D1(:), &
      E0(:)
    real*8, allocatable :: I0(:), Ip(:), IF(:) !< Cumulative integrals of nu_psi parts (0..np)
  end type t_ball_surface

  !> One ballooning problem: either the analytic s-alpha model or a flux surface
  type t_ball_case
    logical :: analytic = .false.
    real*8  :: s_an = 0.d0                 !< Shear for the analytic s-alpha model
    type(t_ball_surface), pointer :: surf => null()
    real*8  :: shear_factor = 1.d0         !< Target shear / actual shear for the surface variation
    logical :: drive_only = .false.        !< Scale P0' only in the curvature drive (no change of
                                           !< the local shear; diagnostic, not self-consistent)
    integer :: np = 512                    !< Points per 2 pi
    integer :: nturns = 5                  !< theta extends over theta0 +- 2 pi nturns
    integer :: ntheta0 = 1                 !< Number of ballooning angles in [0, 2 pi)
  end type t_ball_case

contains



  !> Free the arrays of a surface
  subroutine ball_cleanup_surface(surf)
    type(t_ball_surface), intent(inout) :: surf
    if ( allocated(surf%theta) ) deallocate(surf%theta, surf%R, surf%Z, surf%Jac, surf%BB2,         &
      surf%gpsi2, surf%rho, surf%D0, surf%D1, surf%E0, surf%I0, surf%Ip, surf%IF)
    surf%valid = .false.
    surf%np    = 0
  end subroutine ball_cleanup_surface



  !> P0' that corresponds to a given alpha_miller on this surface
  real*8 function ball_alpha_to_P0p(surf, alpha)
    type(t_ball_surface), intent(in) :: surf
    real*8,               intent(in) :: alpha
    ! alpha = -2 / (2 pi)^2 dV/dpsi r_minor P0'  (as alpha_miller in mod_expression)
    ball_alpha_to_P0p = - alpha * (2.d0*PI)**2 / ( 2.d0 * surf%dV_dpsi * surf%r_minor )
  end function ball_alpha_to_P0p



  !> Determine the local quantities on nsurf flux surfaces from the n=0 part of the current state.
  subroutine ball_surfaces_from_state(node_list, element_list, eq, nsurf, psin_min, psin_max, np,  &
    nsmall, nmaxstep, deltaphi, surfs, ierr)

    type(type_node_list),              intent(in)  :: node_list
    type(type_element_list),           intent(in)  :: element_list
    type(t_equil_state),               intent(in)  :: eq
    integer,                           intent(in)  :: nsurf, np, nsmall, nmaxstep
    real*8,                            intent(in)  :: psin_min, psin_max, deltaphi
    type(t_ball_surface), allocatable, target, intent(out) :: surfs(:)
    integer,                           intent(out) :: ierr

    integer, parameter :: N_EX = 14
    character(len=12), parameter :: names(N_EX) = (/ 'R           ', 'Z           ', 'Psi_N       ', &
      'BR          ', 'BZ          ', 'Btor        ', 'BR_R        ', 'BR_Z        ', 'BZ_R        ', &
      'BZ_Z        ', 'Btor_R      ', 'Btor_Z      ', 'p_prime_loc ', 'rho         ' /)
    type(t_pol_pos_list) :: pol_pos_list
    type(t_tor_pos_list) :: tor_pos_list
    type(t_expr_list)    :: expr_list
    real*8, allocatable  :: result(:,:,:,:), s0(:), sp(:), sF(:)
    real*8  :: R, BR, BZ, Bt, BR_R, BR_Z, BZ_R, BZ_Z, Bt_R, Bt_Z, psi_R, psi_Z, psi_RR, psi_RZ,       &
      psi_ZZ, gpsi, nR, nZ, tR, tZ, F, F_R, F_Z, Bdott, nu_l, qabs, S, nHn, tHt, h, KR, Kphi, KZ,    &
      a0(3), Bv(3), gv(3), Kv(3), axB(3), gxB(3), dq_over_F_dpsi, dphi
    integer :: i, j, jp, jm, iex, n_flip, k, nmaxs
    real*8  :: ffp(0:np-1)
    type(t_ball_surface), pointer :: bs

    ierr = 0

    ! --- determine_theta_mag interpolates the traced field line linearly onto the theta* grid, so the
    ! --- stored points must not be coarser than the grid, and the number of steps must cover a
    ! --- poloidal turn for q up to 10. (trace_fieldlines restarts each small step from the stored
    ! --- point, so consecutive stored points are deltaphi/nsmallsteps apart in phi.)
    dphi  = min(deltaphi, 2.d0 * PI * real(nsmall) / real(np))
    nmaxs = max(nmaxstep, nint(2.d0 * PI * 10.d0 * real(nsmall) / dphi))
    pol_pos_list = pol_pos(node_list, element_list, eq, nPsiN=nsurf, nTht=np, PsiNmin=psin_min,    &
      PsiNmax=psin_max, nsmallsteps=nsmall, nmaxsteps=nmaxs, deltaphi=dphi)
    if ( .not. allocated(pol_pos_list%pos) ) then
      ierr = 1
      return
    end if
    tor_pos_list = tor_pos(phi=0.d0)

    ! --- Only the n=0 harmonic: all derivatives below are those of the axisymmetric field
    expr_list = exprs(names, N_EX)
    if ( expr_list%n_expr /= N_EX ) then
      write(*,*) 'ERROR in '//trim(THIS_MOD_NAME)//': required expressions missing.'
      ierr = 2
      return
    end if
    call eval_expr(eq, JOREK_UNITS, expr_list, pol_pos_list, tor_pos_list, result, ierr,            &
      only_n0=.true.)
    if ( ierr /= 0 ) return

    allocate( surfs(nsurf), s0(0:np), sp(0:np), sF(0:np) )
    h = 2.d0 * PI / real(np)

    do i = 1, nsurf
      bs => surfs(i)
      bs%np = np
      allocate( bs%theta(0:np-1), bs%R(0:np-1), bs%Z(0:np-1), bs%Jac(0:np-1), bs%BB2(0:np-1),       &
        bs%gpsi2(0:np-1), bs%rho(0:np-1), bs%D0(0:np-1), bs%D1(0:np-1), bs%E0(0:np-1),              &
        bs%I0(0:np), bs%Ip(0:np), bs%IF(0:np) )

      bs%valid   = .not. any( pol_pos_list%pos(1:np,i)%outside )
      bs%volume  = pol_pos_list%pos(1,i)%volume
      bs%dV_dpsi = pol_pos_list%pos(1,i)%dV_dpsi
      bs%r_minor = pol_pos_list%pos(1,i)%r_minor
      dq_over_F_dpsi = pol_pos_list%pos(1,i)%dq_over_F_dpsi
      bs%psin    = sum( result(1,1:np,i,3) ) / real(np)
      bs%P0p     = sum( result(1,1:np,i,13) ) / real(np)
      bs%P0p_var = sqrt( sum( (result(1,1:np,i,13) - bs%P0p)**2 ) / real(np) )
      if ( (.not. bs%valid) .or. (bs%dV_dpsi == 0.d0) ) then
        bs%valid = .false.
        cycle
      end if

      n_flip = 0
      do j = 0, np-1
        k  = j + 1                      ! index in the position list (theta* = 2 pi j / np)
        jp = modulo(j+1, np) + 1
        jm = modulo(j-1, np) + 1
        R    = result(1,k,i,1)
        BR   = result(1,k,i,4);  BZ   = result(1,k,i,5);  Bt   = result(1,k,i,6)
        BR_R = result(1,k,i,7);  BR_Z = result(1,k,i,8)
        BZ_R = result(1,k,i,9);  BZ_Z = result(1,k,i,10)
        Bt_R = result(1,k,i,11); Bt_Z = result(1,k,i,12)
        bs%theta(j) = h * real(j)
        bs%R(j)     = R
        bs%Z(j)     = result(1,k,i,2)
        bs%rho(j)   = result(1,k,i,14)

        ! --- Poloidal flux and its Hessian from B (BR = psi_Z/R, BZ = -psi_R/R, valid for n=0)
        psi_R  = - R * BZ
        psi_Z  =   R * BR
        psi_RR = - BZ - R * BZ_R
        psi_ZZ =   R * BR_Z
        psi_RZ = 0.5d0 * ( - R * BZ_Z + BR + R * BR_R )
        ffp(j) = - ( psi_RR - psi_R / R + psi_ZZ ) - R**2 * bs%P0p
        gpsi   = sqrt( psi_R**2 + psi_Z**2 )
        nR     = psi_R / gpsi
        nZ     = psi_Z / gpsi
        ! --- Unit tangent in the direction of increasing theta*
        tR     = - nZ
        tZ     =   nR
        if ( tR * (pol_pos_list%pos(jp,i)%R - pol_pos_list%pos(jm,i)%R)                             &
           + tZ * (pol_pos_list%pos(jp,i)%Z - pol_pos_list%pos(jm,i)%Z) < 0.d0 ) then
          tR = - tR
          tZ = - tZ
          if ( j > 0 ) n_flip = n_flip + 1
        end if
        nHn = nR*nR*psi_RR + 2.d0*nR*nZ*psi_RZ + nZ*nZ*psi_ZZ
        tHt = tR*tR*psi_RR + 2.d0*tR*tZ*psi_RZ + tZ*tZ*psi_ZZ

        ! --- Toroidal field function and field line pitch
        F     = R * Bt
        F_R   = Bt + R * Bt_R
        F_Z   = R * Bt_Z
        Bdott = BR * tR + BZ * tZ
        nu_l  = ( Bt / R ) / Bdott
        qabs  = abs( F ) * pol_pos_list%pos(1,i)%q_over_F
        if ( j == 0 ) bs%q = sign( qabs, nu_l )

        bs%gpsi2(j) = gpsi**2
        bs%BB2(j)   = BR**2 + BZ**2 + Bt**2
        bs%Jac(j)   = qabs * R**2 / abs(F)

        ! --- Integrands for nu_psi (d theta*): actual state, P0' variation, FF' variation
        S     = ( nR*F_R + nZ*F_Z ) / F - nR / R - nHn / gpsi + tHt / gpsi
        s0(j) = bs%q / gpsi * S
        sp(j) = bs%q / gpsi * R**2 / gpsi
        sF(j) = bs%q / gpsi * ( gpsi / F**2 + 1.d0 / gpsi )

        ! --- Vectors in the right-handed basis (R, phi, Z)
        Bv  = (/ BR, Bt, BZ /)
        gv  = (/ psi_R, 0.d0, psi_Z /)
        a0  = (/ - nu_l * tR, 1.d0 / R, - nu_l * tZ /)  ! grad(alpha) without the nu_psi part
        ! --- (B.grad)B = B^2 kappa + parallel part (which drops out in the products below)
        KR   = BR * BR_R + BZ * BR_Z - Bt**2 / R
        KZ   = BR * BZ_R + BZ * BZ_Z
        Kphi = BR * Bt_R + BZ * Bt_Z + Bt * BR / R
        Kv   = (/ KR, Kphi, KZ /)
        axB  = cross(a0, Bv)
        gxB  = cross(gv, Bv)
        bs%D0(j) = dot_product(axB, Kv)
        bs%D1(j) = dot_product(gxB, Kv)
        bs%E0(j) = dot_product(axB, gv)
      end do
      if ( n_flip > 0 ) then
        write(*,'(a,f8.4,a,i5)') ' WARNING (ballooning): tangent orientation changes, Psi_N=',     &
          bs%psin, ', flips: ', n_flip
        bs%valid = .false.
      end if

      ! --- Cumulative trapezoidal integrals over one period (spectrally accurate at theta = 2 pi)
      bs%FFp     = sum( ffp ) / real(np)
      bs%FFp_var = sqrt( sum( (ffp - bs%FFp)**2 ) / real(np) )                                    &
                 / ( abs(bs%FFp) + sum( bs%R**2 ) / real(np) * abs(bs%P0p) )

      s0(np) = s0(0);  sp(np) = sp(0);  sF(np) = sF(0)
      bs%I0(0) = 0.d0; bs%Ip(0) = 0.d0; bs%IF(0) = 0.d0
      do j = 1, np
        bs%I0(j) = bs%I0(j-1) + 0.5d0 * h * ( s0(j-1) + s0(j) )
        bs%Ip(j) = bs%Ip(j-1) + 0.5d0 * h * ( sp(j-1) + sp(j) )
        bs%IF(j) = bs%IF(j-1) + 0.5d0 * h * ( sF(j-1) + sF(j) )
      end do

      bs%dq_dpsi    = bs%I0(np) / (2.d0 * PI)
      ! --- Finite difference estimate (F constant across surfaces assumed, as for alpha_crit)
      bs%dq_dpsi_fd = sign(1.d0, bs%q) * abs(F) * dq_over_F_dpsi
      bs%alpha0     = - 2.d0 / (2.d0 * PI)**2 * bs%dV_dpsi * bs%r_minor * bs%P0p
    end do

    deallocate( result, s0, sp, sF )
    call cleanup_pol_pos(pol_pos_list)
    call cleanup_tor_pos(tor_pos_list)

  end subroutine ball_surfaces_from_state



  !> Cross product
  pure function cross(a, b) result(c)
    real*8, intent(in) :: a(3), b(3)
    real*8 :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross



  !> Coefficients g, c, w on the extended grid theta = theta0 + h i, i = -M..M
  !! For a surface, alpha is the target alpha_miller (P0' scaled self-consistently, see header).
  subroutine build_coefficients(cs, alpha, k0, M, g, c, w, nu_out)
    type(t_ball_case), intent(in)  :: cs
    real*8,            intent(in)  :: alpha
    integer,           intent(in)  :: k0, M
    real*8,            intent(out) :: g(-M:M), c(-M:M), w(-M:M)
    real*8, optional,  intent(out) :: nu_out(-M:M)

    real*8  :: h, th, th0, lam, P0p_t, dP, dFF, C0, Ct, nu, nu0, ga2
    integer :: i, k, p, turn, np

    np = cs%np
    h  = 2.d0 * PI / real(np)

    if ( cs%analytic ) then
      ! --- s-alpha model, lecture 26 Eq. 26.86
      th0 = h * real(k0)
      do i = -M, M
        th   = th0 + h * real(i)
        lam  = cs%s_an * (th - th0) - alpha * ( sin(th) - sin(th0) )
        g(i) = 1.d0 + lam**2
        w(i) = g(i)
        c(i) = alpha * ( cos(th) + lam * sin(th) )
        if ( present(nu_out) ) nu_out(i) = lam
      end do
      return
    end if

    associate( bs => cs%surf )
      P0p_t = ball_alpha_to_P0p(bs, alpha)
      dP    = P0p_t - bs%P0p
      if ( cs%drive_only ) dP = 0.d0        ! P0p_t still enters the drive c below
      C0    = bs%I0(np)
      Ct    = C0 * cs%shear_factor
      dFF   = ( Ct - C0 - dP * bs%Ip(np) ) / bs%IF(np)
      nu0   = bs%I0(k0) + dP * bs%Ip(k0) + dFF * bs%IF(k0)
      do i = -M, M
        k    = k0 + i
        p    = modulo(k, np)
        turn = (k - p) / np
        nu   = bs%I0(p) + dP * bs%Ip(p) + dFF * bs%IF(p) + real(turn) * Ct - nu0
        ga2  = bs%BB2(p) / bs%gpsi2(p) + nu**2 * bs%gpsi2(p)
        g(i) = ga2 / ( bs%Jac(p) * bs%BB2(p) )
        w(i) = bs%Jac(p) * bs%rho(p) * ga2 / bs%BB2(p)
        c(i) = 2.d0 * bs%Jac(p) * ( bs%D0(p) - nu * bs%D1(p) ) * bs%E0(p) * P0p_t / bs%BB2(p)**3
        if ( present(nu_out) ) nu_out(i) = nu
      end do
    end associate

  end subroutine build_coefficients



  !> Symmetrised tridiagonal matrix of the interior points i = -M+1..M-1 (X(+-M) = 0)
  subroutine tridiagonal(M, h, g, c, w, d, e)
    integer, intent(in)  :: M
    real*8,  intent(in)  :: h, g(-M:M), c(-M:M), w(-M:M)
    real*8,  intent(out) :: d(2*M-1), e(2*M-1)
    real*8  :: gp, gm
    integer :: i, k
    do i = -M+1, M-1
      k  = i + M
      gp = 0.5d0 * ( g(i) + g(i+1) ) / h**2
      gm = 0.5d0 * ( g(i) + g(i-1) ) / h**2
      d(k) = ( - gp - gm + c(i) ) / w(i)
      e(k) = gp / sqrt( w(i) * w(i+1) )   ! e(k) couples k and k+1; e(2M-1) unused
    end do
  end subroutine tridiagonal



  !> Number of eigenvalues of the symmetric tridiagonal matrix (d, e) larger than mu
  integer function count_above(n, d, e, mu)
    integer, intent(in) :: n
    real*8,  intent(in) :: d(n), e(n), mu
    real*8  :: piv
    integer :: k, n_below
    n_below = 0
    piv = d(1) - mu
    if ( piv < 0.d0 ) n_below = 1
    do k = 2, n
      if ( abs(piv) < 1.d-280 ) piv = -1.d-280
      piv = d(k) - mu - e(k-1)**2 / piv
      if ( piv < 0.d0 ) n_below = n_below + 1
    end do
    count_above = n - n_below
  end function count_above



  !> Largest eigenvalue by bisection with Sturm counts
  real*8 function lambda_max_tri(n, d, e)
    integer, intent(in) :: n
    real*8,  intent(in) :: d(n), e(n)
    real*8  :: lo, hi, mid, ee
    integer :: k, it
    hi = -huge(1.d0)
    lo =  huge(1.d0)
    do k = 1, n
      ee = 0.d0
      if ( k > 1 ) ee = ee + abs(e(k-1))
      if ( k < n ) ee = ee + abs(e(k))
      hi = max(hi, d(k) + ee)
      lo = min(lo, d(k) - ee)
    end do
    do it = 1, 200
      mid = 0.5d0 * ( lo + hi )
      if ( count_above(n, d, e, mid) > 0 ) then
        lo = mid
      else
        hi = mid
      end if
      if ( hi - lo <= 1.d-13 * max(abs(lo), abs(hi), 1.d-300) ) exit
    end do
    lambda_max_tri = 0.5d0 * ( lo + hi )
  end function lambda_max_tri



  !> Stability of one case at a given alpha, maximised over the ballooning angles theta0.
  !! n_unst: number of unstable eigenvalues (0 = stable), lam: largest eigenvalue (only computed
  !! if want_lambda), k0_worst: theta0 index of the most unstable angle.
  subroutine ball_evaluate(cs, alpha, n_unst, lam, k0_worst, want_lambda)
    type(t_ball_case), intent(in)  :: cs
    real*8,            intent(in)  :: alpha
    integer,           intent(out) :: n_unst, k0_worst
    real*8,            intent(out) :: lam
    logical,           intent(in)  :: want_lambda

    real*8, allocatable :: g(:), c(:), w(:), d(:), e(:)
    real*8  :: h, lam_j
    integer :: M, n, j, k0, nu_j

    M = cs%nturns * cs%np
    n = 2*M - 1
    h = 2.d0 * PI / real(cs%np)
    allocate( g(-M:M), c(-M:M), w(-M:M), d(n), e(n) )

    n_unst   = 0
    lam      = -huge(1.d0)
    k0_worst = 0
    do j = 0, cs%ntheta0 - 1
      k0 = modulo( nint( real(j) * real(cs%np) / real(cs%ntheta0) ), cs%np )
      call build_coefficients(cs, alpha, k0, M, g, c, w)
      call tridiagonal(M, h, g, c, w, d, e)



      nu_j = count_above(n, d, e, 0.d0)
      if ( want_lambda ) then
        lam_j = lambda_max_tri(n, d, e)
        if ( lam_j > lam ) then
          lam      = lam_j
          k0_worst = k0
        end if
      else if ( nu_j > n_unst ) then
        k0_worst = k0
      end if
      n_unst = max(n_unst, nu_j)
    end do

    deallocate( g, c, w, d, e )

  end subroutine ball_evaluate



  !> Marginal alpha values: scan alpha in (0, alpha_max] with nscan steps and bisect the
  !! first unstable (alpha_c1) and the following stable (alpha_c2, second stability) transition.
  !! A value of -1 means that no transition was found in the scanned range.
  subroutine ball_alpha_crit(cs, alpha_max, nscan, alpha_c1, alpha_c2)
    type(t_ball_case), intent(in)  :: cs
    real*8,            intent(in)  :: alpha_max
    integer,           intent(in)  :: nscan
    real*8,            intent(out) :: alpha_c1, alpha_c2

    real*8  :: a_prev, a
    logical :: unst_prev, unst
    integer :: k

    alpha_c1  = -1.d0
    alpha_c2  = -1.d0
    a_prev    = 0.d0
    unst_prev = .false.
    do k = 1, nscan
      a    = alpha_max * real(k) / real(nscan)
      unst = is_unstable(a)
      if ( unst .and. (.not. unst_prev) .and. (alpha_c1 < 0.d0) ) then
        alpha_c1 = bisect(a_prev, a, .false.)
      else if ( (.not. unst) .and. unst_prev .and. (alpha_c1 >= 0.d0) ) then
        alpha_c2 = bisect(a_prev, a, .true.)
        exit
      end if
      a_prev    = a
      unst_prev = unst
    end do

  contains

    logical function is_unstable(alpha)
      real*8, intent(in) :: alpha
      integer :: n_unst, k0
      real*8  :: lam
      call ball_evaluate(cs, alpha, n_unst, lam, k0, .false.)
      is_unstable = ( n_unst > 0 )
    end function is_unstable

    !> Transition between lo (state unst_lo) and hi
    real*8 function bisect(lo_in, hi_in, unst_lo)
      real*8,  intent(in) :: lo_in, hi_in
      logical, intent(in) :: unst_lo
      real*8  :: lo, hi, mid
      integer :: it
      lo = lo_in
      hi = hi_in
      do it = 1, 40
        mid = 0.5d0 * ( lo + hi )
        if ( is_unstable(mid) .eqv. unst_lo ) then
          lo = mid
        else
          hi = mid
        end if
        if ( hi - lo < 1.d-6 * max(hi, 1.d-12) ) exit
      end do
      bisect = 0.5d0 * ( lo + hi )
    end function bisect

  end subroutine ball_alpha_crit



  !> Eigenfunction of the largest eigenvalue for a given alpha and theta0 index (inverse iteration).
  !! Returns theta, X (max |X| = 1), the coefficients, and nu_psi (Lambda for the s-alpha model).
  subroutine ball_eigenfunction(cs, alpha, k0, lam, theta, X, g, c, w, nu)
    type(t_ball_case),   intent(in)  :: cs
    real*8,              intent(in)  :: alpha
    integer,             intent(in)  :: k0
    real*8,              intent(out) :: lam
    real*8, allocatable, intent(out) :: theta(:), X(:), g(:), c(:), w(:), nu(:)

    real*8, allocatable :: d(:), e(:), u(:), piv(:)
    real*8  :: h, sigma
    integer :: M, n, i, k, it

    M = cs%nturns * cs%np
    n = 2*M - 1
    h = 2.d0 * PI / real(cs%np)
    allocate( theta(-M:M), X(-M:M), g(-M:M), c(-M:M), w(-M:M), nu(-M:M), d(n), e(n), u(n), piv(n) )

    call build_coefficients(cs, alpha, k0, M, g, c, w, nu)
    call tridiagonal(M, h, g, c, w, d, e)
    lam   = lambda_max_tri(n, d, e)
    sigma = lam + 1.d-8 * max(abs(lam), 1.d-10)  ! T - sigma I is negative definite

    u(:) = 1.d0
    do it = 1, 4
      ! --- Solve (T - sigma I) u_new = u (LDL^T without pivoting is stable here)
      piv(1) = d(1) - sigma
      do k = 2, n
        piv(k) = d(k) - sigma - e(k-1)**2 / piv(k-1)
        u(k)   = u(k) - e(k-1) / piv(k-1) * u(k-1)
      end do
      u(n) = u(n) / piv(n)
      do k = n-1, 1, -1
        u(k) = ( u(k) - e(k) * u(k+1) ) / piv(k)
      end do
      u = u / maxval(abs(u))
    end do

    do i = -M, M
      theta(i) = h * real(k0 + i)
    end do
    X(-M) = 0.d0
    X( M) = 0.d0
    do i = -M+1, M-1
      X(i) = u(i+M) / sqrt( w(i) )
    end do
    X = X / maxval(abs(X))
    if ( X(0) < 0.d0 ) X = -X

    deallocate( d, e, u, piv )

  end subroutine ball_eigenfunction



  !> Largest eigenvalue for target shear s_t and alpha_t (s-alpha diagram of a surface)
  real*8 function ball_salpha_lambda(cs_in, s_actual, s_t, alpha_t)
    type(t_ball_case), intent(in) :: cs_in
    real*8,            intent(in) :: s_actual, s_t, alpha_t
    type(t_ball_case) :: cs
    integer :: n_unst, k0
    real*8  :: lam
    cs = cs_in
    if ( cs%analytic ) then
      cs%s_an = s_t
    else
      cs%shear_factor = s_t / s_actual
    end if
    call ball_evaluate(cs, alpha_t, n_unst, lam, k0, .true.)
    ball_salpha_lambda = lam
  end function ball_salpha_lambda

end module mod_ballooning
