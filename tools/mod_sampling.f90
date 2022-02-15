!> Module to perform different kinds of sampling tricks
!> including sampling from gaussian distributions and 
!> sampling in cylindrical coordinates.
!>
!> Contains functions like CDF and PDF for different distributions,
!> to use with analytical rootfinding methods to sample from these.
!> That has the advantage that we can sample from any distribution for which
!> we know the CDF, albeit at a slow speed. Be careful with the implementation
!> of CDF and PDF functions if there is a hard edge on the domain. It is probably
!> best to generate a continuous extension of that function for negative numbers
!> for instance, to maintain stability of the rootfinding methods.
!>
!> When adding new methods, try to find a reasonable approximation to the function
!> which is easily invertible. This has an enormous impact on the number of
!> iterations needed to converge. Most distributions look like a hyperbolic tangent
!> in their CDF. This can be scaled and distorted without destroying the invertibility
!> of that approximation. See [[sample_chi_squared_3]] for an example.
module mod_sampling
  use mod_rootfinding
  implicit none
  private
  public :: transform_uniform_cylindrical
  public :: boxmueller_transform
  public :: sample_chi_squared_3
  public :: cdf_chi_squared_3 !< for testing purposes
  public :: sample_gaussian
  public :: sample_dist
  public :: thompson_dist
  public :: sample_cosine
  public :: sample_diffuse
  public :: normal_vectors, normal_vectors_frisvad, normal_vectors_naive
  public :: sample_discrete
  public :: sample_piecewise_linear
  public :: cross_product
  public :: sample_uniform_cone
  public :: sample_uniform_sphere
  public :: sample_uniform_sphere_corona_rthetaphi
  public :: sample_uniform_sphere_corona_rcosphi
  !> Switch here which procedure to use by default. The other ones can be found
  !> by their name
  interface normal_vectors
    module procedure normal_vectors_frisvad
  end interface normal_vectors

  interface transform_uniform_cylindrical
    module procedure transform_uniform_cylindrical_scalar
    module procedure transform_uniform_cylindrical_square
  end interface transform_uniform_cylindrical

  interface sample_uniform_cone
    module procedure sample_uniform_standard_cone
    module procedure sample_uniform_direction_cone
    module procedure sample_uniform_direction_length_cone
  end interface sample_uniform_cone

  interface sample_uniform_sphere_corona_rthetaphi
    module procedure sample_uniform_sphere_corona_rthetaphi_r8
  end interface sample_uniform_sphere_corona_rthetaphi

  interface sample_uniform_sphere_corona_rcosphi
    module procedure sample_uniform_sphere_corona_rcosphi_r8
  end interface sample_uniform_sphere_corona_rcosphi

  !> The Thompson distribution with parameters E_b and n
  type, extends(ddfun) :: thompson_dist
    real*8 :: E_b = 10 ! [eV]
    integer :: n = 2 ! high-energy tail fall-of
  contains
    procedure :: inverse_f => guess_inverse_CDF_thompson
    procedure :: f => CDF_thompson
    procedure :: df => PDF_thompson
    procedure :: ddf => DPDF_thompson
  end type

contains
  !> Transform three uniform random numbers in [0,1] to 
  !> Uniform random numbers in cylindrical coordinates (R,Z,Phi)
  pure subroutine transform_uniform_cylindrical_scalar(ran3, Rbox, Zbox, Phibox, R, Z, Phi)
    implicit none
    real*8, dimension(3), intent(in) :: ran3
    real*8, intent(in), dimension(2) :: Rbox, Zbox, Phibox
    real*8, intent(out) :: R, Z, Phi
    ! Use inversion sampling to correct for cylindrical coordinates
    ! r = sqrt(rand() (B^2 - A^2) + A^2) for min and max radius A and B
    R   = sqrt(ran3(1) * (Rbox(2)**2-Rbox(1)**2) + Rbox(1)**2)
    Z   = (Zbox(2)-Zbox(1))*ran3(2) + Zbox(1)
    phi = (Phibox(2)-Phibox(1))*ran3(3) + Phibox(1)
  end subroutine transform_uniform_cylindrical_scalar

  !> Transform three uniform random numbers in [0,1] to 
  !> uniform random numbers in cylindrical coordinates (R,Z,phi)
  !> inputs:
  !>   ran3:   (real8)(3) array of three random numbers within [0,1]
  !>   Rbox2:  (real8)(2) squared of the major radius bounding box
  !>   Zbox:   (real8)(2) bounding box of the vertical position
  !>   phibox: (real8)(2) bounding box of the toroidal angle
  !> outputs:
  !>   RZPhi:  (real8)(3) uniform sample in cilyndrical coordinates
  pure subroutine transform_uniform_cylindrical_square(ran3,Rbox2,Zbox,phibox,RZPhi)
    implicit none
    real*8,dimension(3),intent(in)  :: ran3
    real*8,dimension(2),intent(in)  :: Rbox2,Zbox,phibox
    real*8,dimension(3),intent(out) :: RZPhi
    RZPhi = (/sqrt(Rbox2(1)+(Rbox2(2)-Rbox2(1))*ran3(1)),&
            Zbox(1)+(Zbox(2)-Zbox(1))*ran3(2),&
            phibox(1)+(phibox(2)-phibox(1))*ran3(3)/)
  end subroutine transform_uniform_cylindrical_square

  !> Transform 2N uniform random numbers in [0,1] to
  !> gaussian-distributed random numbers with mean 0 and sigma 1
  !> Using the box-muller method (very slow!)
  pure function boxmueller_transform(ran) result(out)
    implicit none
    real*8, dimension(:), intent(in) :: ran
    real*8, parameter :: TWOPI = 6.2831853071795864769d0 
    real*8, dimension(size(ran,1)) :: out
    integer :: i

    do i=1,size(ran,1),2
      out(i:i+1) = sqrt(-2*log(max(ran(i),1d-200))) * &
          [cos(TWOPI*ran(i+1)), sin(TWOPI*ran(i+1))]
    end do
  end function boxmueller_transform

  !> Transform a uniformly distributed number u on [0,1] into a normally distributed
  !> number by inverse transform sampling (slow!)
  !> It is much better to use box-muller or something else
  pure function sample_gaussian(u) result(x)
    real*8, intent(in)   :: u !< Uniformly distributed number in [0,1]
    real*8               :: x !< Normally distributed number
    real*8               :: x0 !< Initial guess
    integer              :: ierr
    real*8, parameter    :: guess_a = 1.20278251

    ! Generate a guess by inverting a nearby distribution:
    ! u = tanh(a*x)
    x0 = sqrt(2.d0)*(atanh(2.d0*u-1.d0)/guess_a)
    
    call halleys_method(f=CDF_gaussian, &
                       df=PDF_gaussian, &
                      ddf=PDF_prime_gaussian, &
                       y0=u, x0=x0, x=x, ierr=ierr)
    ! Dangerous: ignore ierr for now
  end function sample_gaussian

  !> Sample from a Knudsen cosine distribution, representing well the angle
  !> of sputtered particles
  !> see https://www.sciencedirect.com/science/article/pii/S0042207X02001732
  !> and https://www-sciencedirect-com.dianus.libr.tue.nl/science/article/pii/0168583X89907040
  pure function sample_cosine(u, n) result(v)
    use constants, only: TWOPI
    real*8, intent(in) :: u(2) !< 2 uniformly distributed random numbers
    real*8, intent(in) :: n(3) !< Normal vector of the surface (does not need to be normalized to 1)
    real*8 :: v(3) !< Resulting velocity vector

    real*8 :: sin_theta, cos_theta, psi, a, b, c
    real*8 :: normal(3), tang1(3), tang2(3)

    sin_theta = sqrt(u(1))
    cos_theta = sqrt(1.d0 - u(1))

    psi = TWOPI * u(2)

    a = sin_theta*cos(psi)
    b = sin_theta*sin(psi)
    c = cos_theta

    ! prepare the normal vector
    normal = n/norm2(n)
    ! make vectors perpendicular to the normal vector
    call normal_vectors(normal, tang1, tang2)
    v = c*normal + a*tang1 + b*tang2
  end function sample_cosine

  !> This corresponds to diffuse emission, as the Lambert cosine distribution. See
  !> https://www.particleincell.com/2015/cosine-distribution/
  !> and https://en.wikipedia.org/wiki/Lambert%27s_cosine_law
  !> indicating uniform emission in solid angle
  pure function sample_diffuse(u, n) result(v)
    use constants, only: TWOPI
    real*8, intent(in) :: u(2) !< 2 uniformly distributed random numbers
    real*8, intent(in) :: n(3) !< Normal vector of the surface (does not need to be normalized to 1)
    real*8 :: v(3) !< Resulting velocity vector

    real*8 :: sin_theta, cos_theta, psi, a, b, c
    real*8 :: normal(3), tang1(3), tang2(3)

    sin_theta = u(1)
    cos_theta = sqrt(1.d0 - sin_theta**2)

    psi = TWOPI * u(2)

    a = sin_theta*cos(psi)
    b = sin_theta*sin(psi)
    c = cos_theta

    ! prepare the normal vector
    normal = n/norm2(n)
    ! make vectors perpendicular to the normal vector
    call normal_vectors(normal, tang1, tang2)
    v = c*normal + a*tang1 + b*tang2
  end function sample_diffuse

  !> Construct normal vectors with a naive algorithm
  pure subroutine normal_vectors_naive(a, b, c)
    real*8, dimension(3), intent(in) :: a
    real*8, dimension(3), intent(out) :: b, c
    if (abs(a(1)) > 0.9d0) then
      b = [0.d0, 1.d0, 0.d0]
    else
      b = [1.d0, 0.d0, 0.d0]
    end if
    b = b - a*dot_product(b, a)
    b = b / norm2(b)
    c = cross_product(a, b)
  end subroutine normal_vectors_naive

  !> From http://orbit.dtu.dk/files/126824972/onb_frisvad_jgt2012_v2.pdf
  pure subroutine normal_vectors_frisvad(n, b1, b2)
    real*8, dimension(3), intent(in) :: n !< Unit vector input (must be normalized)
    real*8, dimension(3), intent(out) :: b1, b2
    real*8 :: a, b, a_norm(3)

    if (n(3) < -.9999999d0) then ! Handle the singularity
      b1 = [0.d0, -1.d0, 0.d0]
      b2 = [-1.d0, 0.d0, 0.d0]
      return
    end if
    a = 1.d0/(1.d0 + n(3))
    b = -n(1)*n(2)*a
    b1 = [1.d0 - n(1)*n(1)*a, b, -n(1)]
    b2 = [b, 1.d0 - n(2)*n(2)*a, -n(2)]
  end subroutine normal_vectors_frisvad

  !> Transform a uniformly distributed number u on [0,1] into a chi^2(3)-distributed
  !> number by inverse transform sampling.
  !> This function might be a bit expensive. There are cheaper ways of making
  !> chi-squared distributions from squares of normally distributed variables.
  !> If you are using quasi-random sequences those are not suitable however.
  pure function sample_chi_squared_3(u) result(x)
    real*8, intent(in)   :: u !< Uniformly distributed number in [0,1]
    real*8               :: x !< Chi-squared(3) distributed number
    real*8               :: x0 !< Initial guess
    integer              :: ierr
    real*8, parameter    :: guess_a = 0.24851051d0, guess_b = 1.11289237d0

    ! Generate a guess by inverting a nearby distribution:
    ! u = tanh(a*x)**b
    x0 = atanh(u**(1.d0/guess_b))/guess_a
    
    call newtons_method(f=CDF_chi_squared_3, &
                       df=PDF_chi_squared_3, &
                       y0=u, x0=x0, x=x, ierr=ierr)
    ! Dangerous: ignore ierr for now
  end function sample_chi_squared_3

  pure function PDF_chi_squared_3(x) result(P)
    real*8, intent(in) :: x
    real*8             :: P
    real*8, parameter  :: PI = atan2(0.d0, -1.d0)
    ! Continue for negative values
    P = sign(sqrt(abs(x))*exp(-abs(x)*0.5d0)/sqrt(2.d0*PI), x)
  end function PDF_chi_squared_3

  pure function CDF_chi_squared_3(x) result(C)
    real*8, intent(in) :: x
    real*8             :: C
    real*8, parameter  :: PI = atan2(0.d0, -1.d0)
    ! Continue for negative values, mirror about y=-x
    C = sign(erf(sqrt(abs(x))/sqrt(2.d0)) - sqrt(2.d0/PI)*exp(-abs(x)*0.5d0)*sqrt(abs(x)), x)
  end function CDF_chi_squared_3

  !> Second derivative of the CDF for a chi_squared(3) distribution.
  !> Use this if you would like to use Halley's method to sample.
  !> In my tests it is faster to use Newton-Raphson iteration though.
  pure function PDF_prime_chi_squared_3(x) result(P)
    real*8, intent(in) :: x
    real*8             :: P
    real*8, parameter  :: PI = atan2(0.d0, -1.d0)
    ! Workaround if we go out of domain
    P = sign((1-abs(x))*exp(-abs(x)*0.5d0)/(2.d0*sqrt(2.d0*PI*abs(x))), x)
  end function PDF_prime_chi_squared_3

  !> PDF of a gaussian (normal) distribution with sigma = 1 and mu = 0
  pure function PDF_gaussian(x) result(P)
    use constants, only: PI
    real*8, intent(in) :: x
    real*8             :: P
    P = 1.d0/sqrt(2.d0*PI) * exp(-x**2/2.d0)
  end function PDF_gaussian

  !> Second derivative of a CDF of a gaussian (normal) distribution with sigma = 1 and mu = 0
  pure function PDF_prime_gaussian(x) result(P)
    use constants, only: PI
    real*8, intent(in) :: x
    real*8             :: P
    P = 1.d0/sqrt(2.d0*PI) * exp(-x**2/2.d0) * (-x)
  end function PDF_prime_gaussian

  !> CDF of a gaussian (normal) distribution with sigma = 1 and mu = 0
  pure function CDF_gaussian(x) result(P)
    real*8, intent(in) :: x
    real*8             :: P
    P = 0.5d0*(1.d0 + erf(x/sqrt(2.d0)))
  end function CDF_gaussian

  !> ---- Thompson distribution sampling functions (CDF, PDF, PDF')
  !> \( \mathrm{CDF}(x) = \left[ 1 - E_b^{n-1} \frac{E_b + n x}{(E_b + x)^n} \right] \)
  pure function CDF_thompson(this, x)
    class(thompson_dist), intent(in) :: this
    real*8, intent(in) :: x
    real*8 :: CDF_thompson
    ! Force x to positive and use the sign later to make the continuation to negative numbers
    ! which is needed for stability of the iterative method
    CDF_thompson = 1.d0 - this%E_b**(this%n-1) * (this%E_b + real(this%n,8) * abs(x))/((this%E_b + abs(x))**this%n)
    CDF_thompson = sign(1.d0, x)*CDF_thompson ! if x < 0 inverse the result
  end function CDF_thompson

  !> \(  n(n-1) \frac{E E_b^{n-1}}{{(E + E_b)}^{n+1}} \)
  pure function PDF_thompson(this, x)
    class(thompson_dist), intent(in) :: this
    real*8, intent(in) :: x
    real*8 :: PDF_thompson
    PDF_thompson = real(this%n * (this%n-1), 8) * (abs(x) * this%E_b**(this%n - 1))/((abs(x) + this%E_b)**(this%n+1))
    PDF_thompson = sign(1.d0, x)*PDF_thompson ! if x < 0 inverse the result
  end function PDF_thompson

  !> \(  n(n-1) E_b^{n-1} \frac{E_b - n x}{(E_b + x)^{(n+2)}} \)
  pure function DPDF_thompson(this, x)
    class(thompson_dist), intent(in) :: this
    real*8, intent(in) :: x
    real*8 :: DPDF_thompson
    DPDF_thompson = real(this%n * (this%n-1), 8) * (this%E_b**(this%n - 1)) * &
        (this%E_b - real(this%n,8)*abs(x))/((abs(x) + this%E_b)**(this%n+2))
    DPDF_thompson = sign(1.d0, x)*DPDF_thompson ! if x < 0 inverse the result
  end function DPDF_thompson

  !> We make a very shitty estimate of this function from some scaling arguments
  !> We could do much better by minimizing the integral of some parametrized
  !> invertible function against E_b and n but I have no time for it.
  !> The function we find is f = 1-1/(1+x/E_b/2.3)
  !> we don't take n into account here (n=2 is very common) and 2.3 is really
  !> arbitrary (but fits ok for E_b = 10, 20). This function is accurate to within
  !> 10% or so.
  !>
  !> If n == 2,3 we can solve the n-th order polynomial equation efficiently
  !> so our guess for the inverse is exact.
  pure function guess_inverse_CDF_thompson(this, f) result(x)
    class(thompson_dist), intent(in) :: this
    real*8, intent(in) :: f
    real*8 :: x
    real*8 :: A, B, C, D

    select case (this%n)
    case (2)
      ! root(A,B,C,SGN) for A x^2 + B x + C = 0, D = B^2-4 A C, sgn indicates
      ! which root to take
      ! We write the equation in terms of y = E_b + x
      A = (1.d0 - f)/this%E_b**(this%n-1)
      B = -2.d0
      C = real(this%n-1,8)*this%E_b
      D = B*B-4*A*C

      x = root(A, B, C, D, 1.d0) - this%E_b
    !case (3) ! not implemented yet
    case default
      x = this%E_b * 2.3d0 * (f / (1.d0 - f))
    end select
  end function guess_inverse_CDF_thompson

  !> Transform a uniformly distributed number u on [0,1) into a number distributed
  !> with the passed dist by inverse transform sampling.
  function sample_dist(dist, u) result(x)
    class(fun), intent(in) :: dist
    real*8, intent(in)   :: u !< Uniformly distributed number in [0,1]
    real*8               :: x !< Normally distributed number
    integer              :: ierr

    select type(dist)
    class is (dfun)
      call newtons_method(dist, y0=u, x=x, ierr=ierr)
    class is (ddfun)
      call halleys_method(dist, y0=u, x=x, ierr=ierr)
    end select

    if (ierr .ne. 0) then
      write(*,*) "Error in sampling distribution with u=", u, ierr, " err: ", dist%f(x) - u
      write(*,*) "Guess quality", dist%inverse_f(u), x
    end if
  end function sample_dist

  !> Sample from a discrete distribution, i.e. a list of probabilities.
  !> This is done by calculating the cumulative sum array and bisecting it with our
  !> uniform random number.
  !> See https://stackoverflow.com/questions/16489449/select-element-from-array-with-probability-proportional-to-its-value#16490300
  !> We take a vector of probabilities p and a uniform random number on [0,1] and return
  !> the index into p of the selected value.
  pure function sample_discrete(p, u) result(i_out)
    real*8, dimension(:), intent(in) :: p !< Vector of probabilities of size >1
    real*8, intent(in)               :: u !< Uniformly distributed random number on [0,1]
    integer :: i_out

    integer :: i, i_min, i_max, i_try
    real*8, dimension(:), allocatable :: c !< Cumulative probability density vector
    real*8 :: split !< The value to find

    allocate(c(size(p,1)))
    c(1) = p(1)
    do i=2,size(p,1)
      c(i) = c(i-1) + p(i)
    end do

    split = u*c(size(p,1))
    i_min = lbound(p,1)
    i_max = ubound(p,1)
    do i=1,size(p,1) ! upper bound on number of steps
      i_try = (i_min + i_max)/2
      if (split .gt. c(i_try)) then
        i_min = i_try
      else
        i_max = i_try
      end if
      if (i_max - i_min .le. 1) then
        i_out = i_max
        return
      end if
    end do
    i_out = size(p,1) ! it must be the last one?
  end function sample_discrete

  !> Sample from a linearly interpolated probability density
  !> This is done by calculating the cumulative sum array and bisecting it with our
  !> uniform random number. Then we calculate the root of the quadratic polynomial manually.
  function sample_piecewise_linear(n, x, p, u) result(x_out)
    integer, intent(in)              :: n !< Number of points in domain
    real*8, dimension(n), intent(in) :: x !< Vector of positions of size n
    real*8, dimension(n), intent(in) :: p !< Vector of probability densities at points
    real*8, intent(in)               :: u !< Uniformly distributed random number on [0,1]
    real*8 :: x_out

    integer :: i, i_min, i_max, i_try
    real*8, dimension(:), allocatable :: c !< Cumulative probability density vector
    real*8 :: split !< The value to find
    real*8 :: a, b, d

    allocate(c(n))
    c(1) = 0.d0
    do i=2,size(p,1)
      c(i) = c(i-1) + 0.5d0 * (p(i) + p(i-1)) * (x(i) - x(i-1))
    end do

    split = u*c(ubound(p,1))
    i_min = lbound(p,1)
    i_max = ubound(p,1)
    do i=1,size(p,1) ! upper bound on number of steps
      if (i_max - i_min .eq. 1) then
        a = (p(i_max)-p(i_min))/(x(i_max)-x(i_min))/2
        b = p(i_min)
        d = max(b*b - 4*a*(c(i_min)-split), 0.d0) ! force to zero at edges against numerical issues

        if (abs(a) .gt. 1d-10) then
          x_out = x(i_min) + (-b + sqrt(d))/(2*a) ! always take this root
        else
          x_out = x(i_min) + (split - c(i_min))/p(i_min)
        end if
        return
      end if

      i_try = (i_min + i_max)/2
      if (split .gt. c(i_try)) then
        i_min = i_try
      else
        i_max = i_try
      end if
    end do
  end function sample_piecewise_linear

  !> uniform sampling of the oriented cone modifing the ray length
  !> inputs:
  !>   cos_half_angle: (real8) cosinus of the cone half angle 
  !>   u:              (real8)(2) random numbers
  !>   dir:            (real8)(3) vector providing the direction
  !>   origin:         (real8)(3) new cone origin
  !>   len_int:        (real8)(2) length interval
  !> outputs:
  !>   ray: (real8)(3) point on the directional unit cone
  function sample_uniform_direction_length_cone(cos_half_angle,u,dir,&
  origin,len_int) result(ray)
    implicit none
    real*8,intent(in)              :: cos_half_angle
    real*8,dimension(2),intent(in) :: len_int
    real*8,dimension(3),intent(in) :: u,dir,origin
    real*8,dimension(3)            :: ray
    real*8                         :: length
    real*8,dimension(3)            :: delta
    !> compute oriented cone
    ray = sample_uniform_direction_cone(cos_half_angle,u(1:2),dir,origin)
    !> modify ray length
    length = len_int(1) + (len_int(2)-len_int(1))*u(3)
    delta = ray-origin
    ray = origin + (delta/sqrt(delta(1)*delta(1)+delta(2)*delta(2)+&
    delta(3)*delta(3)))*length
 end function sample_uniform_direction_length_cone

  !> uniform sampling of the oriented cone given its half aperture angle
  !> and a direction. Possible to modify the cone origin has well
  !> inputs:
  !>   cos_half_angle: (real8) cosinus of the cone half angle 
  !>   u:              (real8)(2) random numbers
  !>   dir:            (real8)(3) vector providing the direction
  !>   origin:         (real8)(3) new cone origin
  !> outputs:
  !>   ray: (real8)(3) point on the directional unit cone
  function sample_uniform_direction_cone(cos_half_angle,u,dir,origin) result(ray)
    implicit none
    real*8,intent(in) :: cos_half_angle
    real*8,dimension(2),intent(in) :: u
    real*8,dimension(3),intent(in) :: dir,origin
    real*8,dimension(3)            :: ray
    real*8                         :: cos_alpha,cos_beta !< spherical coord.
    real*8                         :: sin_alpha,sin_beta               
    real*8,dimension(3)            :: t_vec 
    real*8,dimension(3,3)          :: rot_matrix
    !> compute 
    t_vec = dir/sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
    cos_alpha = acos(t_vec(3)); cos_beta = atan2(t_vec(2),t_vec(1));
    sin_alpha = sin(cos_alpha); sin_beta = sin(cos_beta);
    cos_alpha = cos(cos_alpha); cos_beta = cos(cos_beta);
    !> compute random points on the standard sphere
    ray = sample_uniform_standard_cone(cos_half_angle,u)
    !> rotate the array along the dir direction
    rot_matrix(:,1) = (/-sin_beta,cos_beta,0.d0/)
    rot_matrix(:,2) = (/cos_beta*cos_alpha,sin_beta*cos_alpha,-sin_alpha/)
    rot_matrix(:,3) = (/cos_beta*sin_alpha,sin_beta*sin_alpha,cos_alpha/)
    ray = matmul(rot_matrix,ray)
    !> modify cone origin
    ray = origin + ray
  end function sample_uniform_direction_cone 

  !> uniform sampling of the standard cone given its half aperture angle
  !> the cone axis is assumed to be along the unit sphere north pole
  !> the direction is given in cartesian coordinates
  !> inputs:
  !>   cos_half_angle: (real8) cosinus of the half angle
  !>   u:              (real8)(2) random number uniformly distributed
  !> outputs
  !>   ray: (real8)(3) point on the standard unit cone
  function sample_uniform_standard_cone(cos_half_angle,u) result(ray)
    use constants, only: TWOPI
    implicit none
    real*8,intent(in)              :: cos_half_angle
    real*8,dimension(2),intent(in) :: u
    real*8,dimension(3)            :: ray
    real*8                         :: z2
    real*8 ,dimension(2)           :: zphi
    !zphi =(/0.d0,cos_half_angle/) + (/TWOPI,1.d0-cos_half_angle/)*u
    zphi(1) = TWOPI*u(1)
    zphi(2) = cos_half_angle + (1.d0-cos_half_angle)*u(2)
    z2 = sqrt(1.d0-zphi(2)*zphi(2))
    ray = (/z2*cos(zphi(1)),z2*sin(zphi(1)),zphi(2)/)
  end function sample_uniform_standard_cone

  !> uniform sampling in a sphere volume
  !> inputs:
  !>   u:         (real8)(3) uniform random numbers
  !>   R:         (real8) sphere radius
  !>   cos_theta: (real8)(2) interval of the cosinus of the colatitude  in [-1,1]
  !>   phi:       (real8)(2) azimuth interval in [0,2*pi]
  !> outputs: 
  !>   RThetaPhi: (real8)(3) random R,theta,phi coordinates
  function sample_uniform_sphere(R,cos_theta,phi,u) result(RThetaPhi)
    implicit none
    !> inputs:
    real*8,intent(in) :: R
    real*8,dimension(2),intent(in) :: cos_theta,phi
    real*8,dimension(3),intent(in) :: u
    !> outputs:
    real*8,dimension(3) :: RThetaPhi
    !> compute uniform samples
    RThetaPhi = (/R*(u(1)**(1.d0/3.d0)),&
                acos((cos_theta(2)-cos_theta(1))*u(2)+cos_theta(1)),&
                phi(1)+(phi(2)-phi(1))*u(3)/)
  end function sample_uniform_sphere 

  !> sampling sphere corona using angles coordinates
  !> inputs:
  !>   R_cube:    (real8)(2) cube of the sphere radius interval
  !>   cos_theta: (real8)(2) interval of the cosinus of the colatitude in [-1,1]
  !>   phi:       (real8)(2) azimuth interval in [0,2*pi]
  !>   u:         (real8)(3) uniform random numbers in [0,1]
  !> outputs:
  !>   RThetaPhi: (real8)(3) radom radius, colatitude and azimuth
  function sample_uniform_sphere_corona_rthetaphi_r8(R_cube,cos_theta,phi,u) result(RThetaPhi)
    implicit none
    !> inputs:
    real*8,dimension(2),intent(in) :: R_cube,cos_theta,phi
    real*8,dimension(3),intent(in) :: u
    !> outputs:
    real*8,dimension(3) :: RThetaPhi
    !> compute uniform samples
    RThetaPhi = (/(R_cube(1)+(R_cube(2)-R_cube(1))*u(1))**(1.d0/3.d0),&
    acos(cos_theta(1) + (cos_theta(2)-cos_theta(1))*u(2)),&
    phi(1) + (phi(2)-phi(1))*u(3)/)
  end function sample_uniform_sphere_corona_rthetaphi_r8

  !> sampling sphere corona using R, cosinus, azimuth coordinates
  !> inputs:
  !>   R_cube:    (real8)(2) cube of the sphere radius interval
  !>   cos_theta: (real8)(2) interval of the cosinus of the colatitude in [-1,1]
  !>   phi:       (real8)(2) azimuth interval in [0,2*pi]
  !>   u:         (real8)(3) uniform random numbers in [0,1]
  !> outputs:
  !>   RCosPhi: (real8)(3) radom radius, cosinus of the colatitude and azimuth
  function sample_uniform_sphere_corona_rcosphi_r8(R_cube,cos_theta,phi,u) result(RCosPhi)
    implicit none
    !> inputs:
    real*8,dimension(2),intent(in) :: R_cube,cos_theta,phi
    real*8,dimension(3),intent(in) :: u
    !> outputs:
    real*8,dimension(3) :: RCosPhi
    !> compute uniform samples
    RCosPhi = (/(R_cube(1)+(R_cube(2)-R_cube(1))*u(1))**(1.d0/3.d0),&
    cos_theta(1) + (cos_theta(2)-cos_theta(1))*u(2),&
    phi(1) + (phi(2)-phi(1))*u(3)/)
  end function sample_uniform_sphere_corona_rcosphi_r8

  pure function cross_product(a, b)
    real*8, dimension(3) :: cross_product
    real*8, dimension(3), intent(in) :: a, b

    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product
end module mod_sampling
