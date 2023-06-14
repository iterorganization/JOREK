!> Contains an implementation of the binary collision model, which samples a particle
!> from a background plasma described by a temperature, background flow and temperature
!> gradient. See Y. Homma, A. Hatayama / JCP 231 (2012)
module mod_collisions
  implicit none
  private
  public :: collide_particles
  public :: coulomb_logarithm
  public :: q_homma2013
  public :: sample_velocity_dist_magnetized
  public :: sample_velocity_dist_unmagnetized
contains


pure function coulomb_logarithm(kT, n_i, q_a, q_b, m_a, m_b)
  use constants, only: PI, EPS_ZERO, ATOMIC_MASS_UNIT, EL_CHG
  real*8, intent(in) :: kT !< Temperature [J]
  real*8, intent(in) :: n_i !< Dominant ion density [m^-3]
  integer*1, intent(in) :: q_a, q_b !< charge states of trace (a) and background (b) ions
  real*8, intent(in) :: m_a, m_b !< atomic masses of trace and background ions in atomic mass units [u]
  real*8 :: coulomb_logarithm

  coulomb_logarithm = log(sqrt(EPS_ZERO * kT / (n_i * EL_CHG**2)) * &
      4 * PI * EPS_ZERO * ((m_a * m_b)/(m_a + m_b)) * ATOMIC_MASS_UNIT * &
      kT / (m_b * ATOMIC_MASS_UNIT) / &
      (real(q_a,8) * real(q_b,8) * EL_CHG**2))
! coulomb_logarithm = 15.d0 ! for Homma JCP 2013 benchmark
! We should probably use a simpler parametrisation for performance reasons
end function coulomb_logarithm

!> Collide two particles according to the method in Takizuka, Abe 1977
pure subroutine collide_particles(u, q_a, m_a, v_a, q_b, m_b, v_b, n_b, coulomb_logarithm, dt)
  use constants
  use mod_sampling, only: boxmueller_transform

  real*8, dimension(3), intent(in)   :: u !< uniformly distributed random numbers
  integer*1, intent(in)              :: q_a !< charge state [e]
  real*8, intent(in)                 :: m_a !< Ion mass in atomic mass units
  real*8, intent(inout)              :: v_a(3) !< particle velocity [m/s]
  integer*1, intent(in)              :: q_b !< charge state [e]
  real*8, intent(in)                 :: m_b !< Ion mass in atomic mass units
  real*8, intent(inout)              :: v_b(3) !< particle velocity [m/s]
  real*8, intent(in)                 :: n_b !< Background density [m^-3]
  real*8, intent(in)                 :: coulomb_logarithm !< ln(Lambda) indicates how important large or small-angle deflections are. Between 5 and 15 usually
  real*8, intent(in)                 :: dt !< timestep between collisions

  real*8 :: phi !< polar angle phi [radians]
  real*8 :: sin_theta, cos_phi, sin_phi, one_minus_cos_theta  !< goniometric quantities [dimensionless]
  real*8 :: xi(2) !< gaussian random numbers
  real*8 :: Du !< Collisional momentum transfer frequency [1/s]
  real*8 :: v(3), nv, nvp ! relative velocity, norm and perpendicular norm before collision [m/s]
  real*8 :: delta_v(3) ! relative velocity change due to collision [m/s]

  v = v_a - v_b
  nv = norm2(v)
  nvp = norm2(v(1:2))
  Du = real(int(q_a,4)**2 * int(q_b,4)**2,8) * EL_CHG**4 * n_b * coulomb_logarithm / &
      (8.d0 * PI * EPS_ZERO**2 * (m_a*m_b/(m_a+m_b))**2 * ATOMIC_MASS_UNIT**2 * (nv**3))

  xi = boxmueller_transform(u(1:2))*sqrt(Du*dt)
  sin_theta = 2.d0*xi(1)/(1.d0+xi(1)**2)
  one_minus_cos_theta = 2.d0*xi(1)**2/(1.d0+xi(1)**2)
  !theta = 2.d0*atan(xi(1))
  phi   = TWOPI*u(3)
  cos_phi = cos(phi)
  sin_phi = sin(phi)

  if (nvp .lt. 1d-14) then
    delta_v = [nv*sin_theta*cos_phi, nv*sin_theta*sin_phi, -nv*one_minus_cos_theta]
  else
    delta_v = [(v(1)/nvp)*v(3)*sin_theta*cos_phi - (v(2)/nvp)*nv*sin_theta*sin_phi - v(1)*one_minus_cos_theta, &
               (v(2)/nvp)*v(3)*sin_theta*cos_phi + (v(1)/nvp)*nv*sin_theta*sin_phi - v(2)*one_minus_cos_theta, &
               -nvp*sin_theta*cos_phi - v(3)*one_minus_cos_theta]
  end if

  v_a = v_a + (m_b/(m_a+m_b))*delta_v
  v_b = v_b - (m_a/(m_a+m_b))*delta_v
end subroutine collide_particles



!> Calculate the heat flux density q as in the Homma 2013 JCP paper for use with sample_velocity_dist_magnetized
!> The parallel component is in the B_hat direction. The perpendicular component is in the direction
!> of the temperature gradient. The diamagnetic firection is perpendicular to these two.
!>
!> This makes a few assumptions:
!> * Collisional background plasma (no neo or anomalous transport)
!> * Single ion species background plasma
!> * Ion-ion collisions only for background ion
!> * Trace impurity limit
!> * Viscous stress tensor negligible
!> * Strongly magnetized
!> * No pressure gradient or electric field
!>
!> Not all of these assumptions are reasonable everywhere. Be careful.
!>
!> Additionally, it might be better to calculate grad_kT in parallel, dia and perp direction directly
!> instead of reconstructing from real space and the magnetic field vector.
pure function q_homma2013(kT, grad_kT, B, n_b, m_b, q_b) result(q)
  use constants
  use mod_sampling, only: cross_product
  real*8, intent(in) :: kT, grad_kT(3) !< local temperature [J] and gradient in RZPhi space [J/m]
  real*8, intent(in) :: B(3) !< Local magnetic field vector [T]
  real*8, intent(in) :: n_b !< Local background density [m^-3]
  real*8, intent(in) :: m_b !< Background ion mass [u]
  integer*1, intent(in) :: q_b !< Background ion charge state [e]
  real*8 :: q(3) !< Heat flux density vector [W/m^2]

  real*8 :: tau_b !< characteristic collision time of background ion [s]
  real*8 :: K_par, K_dia, K_perp !< parallel, diamagnetic and perpendicular heat conductivities [m^-1 s^-1]
  real*8 :: Omega_b !< Gyrofrequency [rad/s]
  real*8 :: B_hat(3) !< Unit vector in the direction of magnetic field
  real*8 :: grad_kT_perp(3) !< perpendicular gradient of kT [J/m]

  tau_b = 12*pi*sqrt(pi)*EPS_ZERO**2 * sqrt(m_b*ATOMIC_MASS_UNIT) * (kT*sqrt(kT)) / &
      (n_b * real(int(q_b, 4)**4) * EL_CHG**4 * coulomb_logarithm( &
        kT, n_b, q_b, q_b, m_b, m_b &
      ))

  B_hat = B/norm2(B)
  Omega_b = real(q_b,4) * EL_CHG * norm2(B) / (m_b * ATOMIC_MASS_UNIT)

  K_par = 3.9d0 * n_b * kT * tau_b / (m_b * ATOMIC_MASS_UNIT)
  K_dia = 5d0   * n_b * kT / (2.d0 * m_b * ATOMIC_MASS_UNIT * Omega_b)
  K_perp = 2.d0 * n_b * kT / (m_b * ATOMIC_MASS_UNIT * Omega_b**2 * tau_b)

  grad_kT_perp = grad_kT - dot_product(grad_kT, B_hat)*B_hat

  q = (- K_par * dot_product(grad_kT, B_hat) * B_hat &
       + K_dia * cross_product(B_hat, grad_kT_perp) &
       - K_perp * grad_kT_perp)
end function q_homma2013


!> Two variants of this subroutine exist, one for magnetized plasmas and one for unmagnetized plasmas.
!> In the first case we use Homma JCP 2013
!> In the second case we use Homma JCP 2012.
!> Both sampling routines are nearly equivalent. We use a higher-level routine to
!> perform the sampling behind the scenes.
subroutine sample_velocity_dist_magnetized(n, u, ktb, q, n_b, m_b, q_b, v_b, w_1)
  use constants
  integer, intent(in)                :: n !< number of particles to sample
  real*8, dimension(6,n), intent(in) :: u !< uniformly distributed random numbers
  real*8, intent(in)                 :: kTb !< local temperature in [J]
  real*8, intent(in), dimension(3)   :: q !< local heat flux vector in [W/m^2]
  real*8, intent(in)                 :: n_b !< Background density [m^-3]
  real*8, intent(in)                 :: m_b !< Ion mass in atomic mass units
  integer*1, intent(in)              :: q_b !< Background charge state [e]
  real*8, intent(in), dimension(3)   :: v_b !< Background flow velocity [m/s]
  real*8, intent(out), dimension(3,n):: w_1 !< Sampled particle velocities [m/s]
  real*8 :: v_thermal, A
  integer :: i

  v_thermal = sqrt(kTb/(m_b*ATOMIC_MASS_UNIT)) ! [m/s]
  ! A was \(\frac{m_b}{n_b}\frac{1}{(kT_b)^2} \)   [s^4/kg/m]
  ! we multiply A by v_thermal = \(\sqrt{kT_b}{m_b}\) to make w dimensionless
  ! A becomes \(\frac{\sqrt{m_b}}{n_b}\frac{1}{(kT_b)^{3/2}} \) [s^3/kg == m^2/W]
  ! multiplyied by q [W/m^2] this becomes 1
  A = -sqrt(m_b*ATOMIC_MASS_UNIT)/(n_b * kTb*sqrt(kTb)) ! [m^2/W == s^3/kg]
  call sample_distorted_maxwellian(n, u, A, q, w_1)
  ! w_1 is in units of v_thermal now

  ! Add mean flow and normalize with v_thermal
  do i=1,n
    w_1(:,i) = w_1(:,i) * v_thermal + v_b
    if (isnan(w_1(1,i)) .or. isnan(w_1(2,i)) .or. isnan(w_1(3,i))) then
      write(*,*) u(:,i), w_1(:,i)
      call sample_distorted_maxwellian_debug(n, u, A, q, w_1)
    end if
  end do
end subroutine sample_velocity_dist_magnetized



!> sample n particles from the distorted maxwellian with ktb and grad_ktb, and mean flow v_b.
!>
!> If we know m_a and q_a we could also calculate the coulomb_logarithm here instead of presetting it for the simulation
!> (or we do some estimate, an interpolation based on temperature and coronal equilibrium charge state)
!>
!> This routine calculates the heat conductivity of a plasma behind the scenes, and could be replaced
!> by a calculation fo that heat conductivity and [[sample_velocity_dist_magnetized]] which can then
!> be renamed to fit its more general purpose.
subroutine sample_velocity_dist_unmagnetized(n, u, coulomb_logarithm, ktb, grad_ktb, n_b, m_b, q_b, v_b, w_1)
  use constants
  integer, intent(in)                :: n !< number of particles to sample
  real*8, dimension(6,n), intent(in) :: u !< uniformly distributed random numbers
  real*8, intent(in)                 :: coulomb_logarithm !< ln(Lambda) indicates how important large or small-angle deflections are. Between 5 and 15 usually
  real*8, intent(in)                 :: kTb !< local temperature in [J]
  real*8, intent(in), dimension(3)   :: grad_kTb !< local temperature gradient in [J/m]
  real*8, intent(in)                 :: n_b !< Background density [m^-3]
  real*8, intent(in)                 :: m_b !< Ion mass in atomic mass units
  integer*1, intent(in)              :: q_b !< Background charge state [e]
  real*8, intent(in), dimension(3)   :: v_b !< Background flow velocity [m/s]
  real*8, intent(out), dimension(3,n):: w_1 !< Sampled particle velocities, in units of [m/s]
  real*8 :: v_thermal ! [m/s]
  real*8 :: A ! [m/J]
  integer :: i

  v_thermal = sqrt(kTb/(m_b*ATOMIC_MASS_UNIT)) ! [m/s]
  A = (75.d0/2.d0 * PI*sqrt(PI) * EPS_ZERO**2)/coulomb_logarithm * kTb/(n_b*(real(q_b)*EL_CHG)**4) ! [m/J]

  call sample_distorted_maxwellian(n, u, A, grad_kTb, w_1)
  ! Add mean flow and normalize with v_thermal
  do i=1,n
    w_1(:,i) = w_1(:,i) * v_thermal + v_b
  end do
end subroutine sample_velocity_dist_unmagnetized




!> Sample from a distribution of the form
!> 
!> \[ f_b(w) = n_b {\left(\frac{m_b}{2 \pi k T_b}\right)}^{3/2} 
!>      \exp\left( - \frac{m_b w^2}{2 k T_b} \right)
!>      \times \left [1 + A \cdot \left( 1 - \frac{w^2}{5 v_\mathrm{th,b}^2} \right)
!>             \left(q \cdot w\right) \right]
!> \]
!> which is a maxwellian with a skew in a specific direction.
!> 
!> This we simplify by normalizing \(w\) with \(v_\mathrm{th,b}\) and dividing by \(n_b\), yielding
!>
!> \[
!>   f_b(w) = \sqrt{1/2\pi} w^2 \exp\left(-w^2/2\right)
!>        \cdot \left[ 1 + A \left(1-\frac{w^2}{5}\right) (q \cdot w) \right]
!> \]
!> which integrates to 1 exactly for \(q=0\).
!> 
!> We retrieve classical transport by taking
!> \(q = \nabla kT_b\) and \(A = \frac{75}{2} \pi^{3/2} \frac{\epsilon_0^2}{\ln \Lambda} \frac{kT_b}{n_b q_b^4}\). (Homma2012)
!>
!> To obtain neoclassical transport we need a different estimate for the heat flux vector, from neoclassical heat conductivities.
!> The corresponding factor is \(A = - \frac{\sqrt{m_b}}{n_b} {\left(\frac{1}{kT_b}\right)}^{3/2} \) (Homma2013)
!>
!> To automatically calculate these A and q use the routines [[sample_velocity_dist_unmagnetized]] and [[sample_velocity_dist_magnetized]]
pure subroutine sample_distorted_maxwellian(n, u, A, q, w_1)
  use constants, only: TWOPI, EL_CHG
  use mod_sampling

  integer, intent(in)                :: n !< number of particles to sample
  real*8, dimension(6,n), intent(in) :: u !< uniformly distributed random numbers
  real*8, intent(in), dimension(3)   :: q !< local heat flux vector in [W/m^2] (magnetized) or temperature gradient [J/m] (unmagnetized)
  real*8, intent(in)                 :: A !< Distortion coefficient A [m^2/W] (magnetized) or [m/J] (unmagnetized)
  real*8, intent(out), dimension(3,n):: w_1 !< Sampled particle velocities, in units of v_thermal

  real*8 :: w_mag(n), w(4,n), alpha(n), R(n)
  real*8 :: cos_theta_2(n), phi_2(n), sin_theta_2(n)
  real*8 :: T2_1(3,3)
  real*8 :: x, y, z, xy2, rho ! components of grad_kTb, after adding a small x-component
  integer :: i

  ! preparation
  x = sign(max(abs(q(1)), 1d-8*EL_CHG), q(1)) ! add a small component so this works for grad_Ktb = 0 (add 1e-8 eV/m or 1e-8*EL_CHG W/m^2)
  y = q(2)
  z = q(3)
  xy2 = norm2([x,y])
  rho = norm2([x,y,z])
  ! Create transformation matrix in a more efficient way
  T2_1(:,1) = [(x*z)/(xy2*rho), (y*z)/(xy2*rho), -xy2/rho]
  T2_1(:,2) = [-y/xy2, x/xy2, 0.d0]
  T2_1(:,3) = [x/rho,y/rho,z/rho]
  ! Reset the norm for use later
  rho = norm2(q)

  ! 4.2.3
  ! Get 4 uniform random numbers on [0,1]
  ! Transform them into gaussian numbers
  do i=1,n
    w(:,i) = boxmueller_transform(u(1:4,i))
  end do
  ! Get a Chi^2_3-distributed number by summing their squares
  ! and take the square root (here implemented as norm2)
  w_mag = norm2(w(1:3,:),dim=1)! in units of v_thermal

  ! 4.2.4
  ! Removed the divide by v_thermal, multiply later
  alpha = A * (1.d0 - w_mag**2/(5.d0)) * w_mag * rho
  R = u(5,:)
  ! create a single expression for all cases so we don't need if-statements
  ! Expression verified against piecewise expression from paper:
  ! R=0.5, alpha=-1,1
  ! Manipulate[Plot[{Min[(Sqrt[4 Abs[α] R + (1 - Abs[α])^2] - 1)/Abs[α], 2 Sqrt[R] - 1]*Sign[α],
  !  Piecewise[{{2R-1,α==0},{(Sqrt[4 α R + (1-α)^2]-1)/α, 0 < α ≤ 1}, {2Sqrt[R]-1, α>1},
  !                         {(Sqrt[-4α R + (1+α)^2]-1)/α, -1 ≤ α<0}, {1-2Sqrt[R],α< -1}}]
  ! },{α, -1.5, 1.5}], {R,0,1}]
  ! Mistake here with preferential sign of alpha
  cos_theta_2 = 2.d0*R-1.d0 ! default case, where(abs(alpha) .lt. 1d-10)
  where(abs(alpha) .gt. 1.d0) cos_theta_2 = sign(2.d0*sqrt(R)-1.d0,alpha)
  where(abs(alpha) .le. 1.d0 .and. abs(alpha) .gt. 1d-10) cos_theta_2 = (sqrt(4.d0*abs(alpha)*R + (1.d0-abs(alpha))**2) - 1.d0)/alpha

  ! 4.2.5
  phi_2 = TWOPI*u(6,:)
  sin_theta_2 = sqrt(1.d0-min(cos_theta_2**2, 1d0)) ! cap cos_theta_2^2 here
  ! because sometimes it ends up to be something like 1.0000000000002 due to
  ! floating point bullshit

  ! 4.2.6
  do i=1,n
    w_1(:,i) = matmul(T2_1, &
      w_mag(i)*[sin_theta_2(i)*cos(phi_2(i)), sin_theta_2(i)*sin(phi_2(i)), cos_theta_2(i)] &
    )
  end do
end subroutine sample_distorted_maxwellian



subroutine sample_distorted_maxwellian_debug(n, u, A, q, w_1)
  use constants, only: TWOPI, EL_CHG
  use mod_sampling

  integer, intent(in)                :: n !< number of particles to sample
  real*8, dimension(6,n), intent(in) :: u !< uniformly distributed random numbers
  real*8, intent(in), dimension(3)   :: q !< local heat flux vector in [W/m^2]
  real*8, intent(in)                 :: A !< Distortion coefficient A
  real*8, intent(out), dimension(3,n):: w_1 !< Sampled particle velocities, in units of v_thermal

  real*8 :: w_mag(n), w(4,n), alpha(n), R(n)
  real*8 :: cos_theta_2(n), phi_2(n), sin_theta_2(n)
  real*8 :: T2_1(3,3)
  real*8 :: x, y, z, xy2, rho ! components of grad_kTb, after adding a small x-component
  integer :: i

  ! preparation
  x = sign(max(abs(q(1)), 1d-8*EL_CHG), q(1)) ! add a small component so this works for grad_Ktb = 0
  y = q(2)
  z = q(3)
  xy2 = norm2([x,y])
  rho = norm2([x,y,z])
  ! Create transformation matrix in a more efficient way
  T2_1(:,1) = [(x*z)/(xy2*rho), (y*z)/(xy2*rho), -xy2/rho]
  T2_1(:,2) = [-y/xy2, x/xy2, 0.d0]
  T2_1(:,3) = [x/rho,y/rho,z/rho]
  ! Reset the norm for use later
  rho = norm2(q)

  ! 4.2.3
  ! Get 4 uniform random numbers on [0,1]
  ! Transform them into gaussian numbers
  do i=1,n
    w(:,i) = boxmueller_transform(u(1:4,i))
  end do
  ! Get a Chi^2_3-distributed number by summing their squares
  ! and take the square root (here implemented as norm2)
  w_mag = norm2(w(1:3,:),dim=1)! in units of v_thermal

  ! 4.2.4
  ! Removed the divide by v_thermal, multiply later
  alpha = A * (1.d0 - w_mag**2/(5.d0)) * w_mag * rho
  R = u(5,:)
  write(*,*) 'alpha=', alpha
  ! create a single expression for all cases so we don't need if-statements
  ! Expression verified against piecewise expression from paper:
  ! R=0.5, alpha=-1,1
  ! Manipulate[Plot[{Min[(Sqrt[4 Abs[α] R + (1 - Abs[α])^2] - 1)/Abs[α], 2 Sqrt[R] - 1]*Sign[α],
  !  Piecewise[{{2R-1,α==0},{(Sqrt[4 α R + (1-α)^2]-1)/α, 0 < α ≤ 1}, {2Sqrt[R]-1, α>1},
  !                         {(Sqrt[-4α R + (1+α)^2]-1)/α, -1 ≤ α<0}, {1-2Sqrt[R],α< -1}}]
  ! },{α, -1.5, 1.5}], {R,0,1}]
  ! Mistake here with preferential sign of alpha
  cos_theta_2 = 2.d0*R-1.d0 ! default case, where(abs(alpha) .lt. 1d-10)
  where(abs(alpha) .gt. 1.d0) cos_theta_2 = sign(2.d0*sqrt(R)-1.d0,alpha)
  where(abs(alpha) .le. 1.d0 .and. abs(alpha) .gt. 1d-10) cos_theta_2 = (sqrt(4.d0*abs(alpha)*R + (1.d0-abs(alpha))**2) - 1.d0)/alpha
  write(*,*) 'cos(theta)=', cos_theta_2

  ! 4.2.5
  phi_2 = TWOPI*u(6,:)
  write(*,*) min(cos_theta_2**2, 1d0)
  write(*,*) 'sin(theta)^2=', 1.d0-min(cos_theta_2**2, 1d0)
  write(*,*) 'sin(theta)=', sqrt(1.d0-min(cos_theta_2**2, 1d0))
  sin_theta_2 = sqrt(1.d0-min(cos_theta_2**2, 1d0)) ! cap cos_theta_2 here
  ! because sometimes it ends up to be -1.0000000000002 due to floating point
  ! bullshit

  ! 4.2.6
  do i=1,n
    w_1(:,i) = matmul(T2_1, &
      w_mag(i)*[sin_theta_2(i)*cos(phi_2(i)), sin_theta_2(i)*sin(phi_2(i)), cos_theta_2(i)] &
    )
  end do
end subroutine sample_distorted_maxwellian_debug
end module mod_collisions
