!> Mod Eckstein sputter yield and sputtered energy coefficients
!> This module contains types for the needed coefficients for the Eckstein angle
!> dependent fit formulas, for the sputter yield and the sputtered energy coefficient.
!> It evaluates the fit formula at available energies from the list of the eckstein data
!> It interpolates logarithmically between two sputter yields of the closest available energies to the energy of the particle
!> This is done to approximate the real sputter yield belonging to the incoming particle

!> NOTE:
!> The exact same thing is done for the sputtered energy coefficient
!> This is possible as the fitting formula is exactly the same,
!> however, the values for the parameters differ.
!> Instead of the sputter yiled Yn for sputtering from normal incidence,
!> the sputtered energy coefficient for normal incidence is used.
!>
!> See internship report of Sven Korving for more details.
module mod_eckstein_y_ye
use mod_atomic_elements
use constants
 
implicit none
private

public :: eckstein_coeff_set, eckstein_sputter_yield, eckstein_sputtered_energy_coeff
 

!> Eckstein angle dependency fit coefficients, determining the function
!> \[
!>   \frac{Y(E_0,\theta_0)}{Y(E_0,0)} = \left( \cos\left[\left(\frac{\pi \theta}{2 theta^*}\right)^c\right] \right)^{-f}
!>     \exp\left[b\left(1 - \left(\cos\left[\left(\frac{\pi\theta_0}{2\theta^*}\right)^c\right]\right)^{-1}\right]
!> \]
!> where $\theta$ is in degrees, but $\theta^*$ is in radians.
!> See [[evaluate_eckstein_formula]] for an implementation
type :: eckstein_coeff
  real*8 :: E0 !< energy [eV]
  real*8 :: f
  real*8 :: b
  real*8 :: c
  real*8 :: Yn !< Sputtering yield/energy coeff at normal incidence at E0
  real*8 :: Esp !< [eV]
  real*8 :: theta_star !< corrected angle of incidence [degrees]
  real*8 :: theta_max !< angle of maximum sputtering yield [degrees]
  contains
    procedure :: eval => evaluate_eckstein_formula
end type eckstein_coeff



!> Base type for a set of eckstein coefficients for sputter yield or sputtered energy coefficient
type, abstract :: eckstein_coeff_set
  integer :: Z_ion !< atomic number of impacting ion
  integer :: Z_target !< atomic number of target ion

  !> Eckstein coeffs for normal incidence, shared between fits for sputter yield and sputtered energy coefficient
  real*8 :: eps
  real*8 :: q
  real*8 :: mu
  real*8 :: lambda
  real*8 :: E_threshold

  logical :: use_Yn_func = .true. !< By default use the coefficients above to estimate Yn instead of the discrete tabulated data
  
  type(eckstein_coeff), allocatable :: yn(:)
contains
  procedure :: read !< Read yn(:) and other parameters from a file
  procedure(calc_E), deferred, private :: calc_E !< calculate normal incidence yield/coeff
  !< The normal incidence coefficient can be either interpolated from the Yn in the table of yn(:)
  !< or it can be calculated from the coefficients in this set
  procedure, public :: interp !< Call interp_E and interp_theta internally as necessary
end type eckstein_coeff_set

!> Sputter yield calculation
type, extends(eckstein_coeff_set) :: eckstein_sputter_yield
contains
  procedure :: calc_E => calc_E_sputter_yield
end type eckstein_sputter_yield

!> Sputtered energy coefficient calculation
type, extends(eckstein_coeff_set) :: eckstein_sputtered_energy_coeff
  real*8 :: nu
  real*8 :: n
contains
  procedure :: calc_E => calc_E_sputtered_energy_coeff
end type eckstein_sputtered_energy_coeff



interface
  pure function calc_E(this, E) result(yield)
    import eckstein_coeff_set
    class(eckstein_coeff_set), intent(in)  :: this
    real*8, intent(in) :: E
    real*8 :: yield
  end function calc_E
end interface

contains 



!> Interpolate sputter yield or sputtered energy coefficient from components
pure function interp(this, E, theta) result(yield)
  class(eckstein_coeff_set), intent(in) :: this
  real*8, intent(in)                    :: E !< energy in eV
  real*8, intent(in)                    :: theta !< angle in degrees
  real*8 :: yield
  integer :: i1, i2 !< positions to interpolate
  real*8 :: a !< interpolation factor
  if (.not. allocated(this%yn)) then
    yield = 0.d0
    return
  end if

  call interp_factors(E, this%yn, i1, i2, a)

  if (E .le. this%E_threshold) then
    yield = 0
  else
    if (this%use_Yn_func) then
      yield = this%calc_E(E)
    else ! interpolate log-linear
      yield = this%yn(i1)%Yn**a * this%yn(i2)%Yn**(1.d0-a)
    end if
  end if
  ! moving the interpolations is allowed since multiplication is commutative.
  ! i.e. (y_1 theta_1)^a (y_2 theta_2)^(1-a) == y_1^a y_2^(1-a) theta_1^a theta_2^(1-a)
  yield = yield * (this%yn(i1)%eval(theta)**a * this%yn(i2)%eval(theta)**(1.d0-a))
end function interp


!> The common part in the interpolation routines is finding the indices of the
!> reference points and calculating the interpolation factor a. Split these out
!> into this subroutine
pure subroutine interp_factors(E,yn,i1,i2,a)
  real*8, intent(in)   :: E
  type(eckstein_coeff), dimension(:), intent(in) :: yn
  integer, intent(out) :: i1, i2
  real*8, intent(out)  :: a
  integer :: i_delta_E
  real*8 :: E1, E2
  real*8 :: delta_E(size(yn,1))

  delta_E = E - yn(:)%E0
  !> find the two closest (upper and lower) energies that exists in the Eckstein data file for that particle ion target combination
  i_delta_E = minloc(delta_E, dim = 1, mask=delta_E .ge. 0.d0) ! assumes yn(:)%E0 is increasing

  if (i_delta_E .eq. 0) then !< extrapolate downwards
    i1 = 1 !< y1 is the set of coeffs need for evaluate_eckstein_formula at a certain energy
    i2 = 2
  else if (i_delta_E .eq. size(yn,1)) then !< extrapolate upwards
    ! be careful extrapolating in log-space... you could very quickly get large values if you extrapolate too far
    i1 = i_delta_E-1
    i2 = i_delta_E
  else ! default case
    i1 = i_delta_E
    i2 = i_delta_E+1
  end if

  E1 = yn(i1)%E0
  E2 = yn(i2)%E0
  a = (log(E2)-log(E))/(log(E2)-log(E1))
  if (a .lt. 0) then
    a = 0 ! disable extrapolation!
  end if
end subroutine interp_factors


!> Sputter yield fit formula from Eckstein
!> \[
!>   Y(E_0,0) = q S_n^\mathrm{KrC}(\epsilon_L) \frac{\left(\frac{E_0}{E_\mathrm{th}} - 1\right)^\mu}{\lambda + {\left(\frac{E_0}{E_\mathrm{th}} - 1\right)^\mu}}
!> \]
!> at normal angle $\theta = 0$ deg.
!> Here the nuclear stopping power is
!> \[
!>   S_n^\mathrm{KrC}(\epsilon_L) = \frac{0.5 \ln(1+1.2288 \epsilon_L)}{\epsilon_L + 0.1728 \sqrt{\epsilon_L} + 0.008 \epsilon_L^{0.1504}}
!> \]
!> where the reduced energy is
!> \[
!>  \epsilon_L = E_0 \frac{M_2}{M_1 + M_2} \frac{a_L}{Z_1 Z_2 e^2} = \frac{E_0}{\epsilon}
!> \]
!> and \(\epsilon\) is precalculated. See Korving 2018 for details.
pure function calc_E_sputter_yield(this, E) result(yield)
  class(eckstein_sputter_yield), intent(in) :: this
  real*8, intent(in)                        :: E !< energy in eV
  real*8 :: yield
  real*8 :: f_stopping_power !< The nuclear stopping power
  real*8 :: eps_l !< Lindhard reduced energy
  real*8 :: w
  if (.not. allocated(this%yn)) then
    yield = 0.d0
    return
  end if

  eps_l = E/this%eps
  w = eps_l + 0.1728d0 * sqrt(eps_l) + 0.008d0 * (eps_l)**0.1504d0
  
  f_stopping_power = 0.5d0 * log(1.d0+1.2288d0*eps_l) / w

  yield = this%q * f_stopping_power * (E/this%E_threshold - 1.d0)**this%mu / (this%lambda/w + (E/this%E_threshold - 1.d0)**this%mu)
end function calc_E_sputter_yield

!> Sputtered energy coefficient fit formula. Same as the sputter yield but with extra power $n$ and $\nu$ in stopping power.
!> \[
!>   S_n^\mathrm{KrC}(\epsilon_L) = \left(\frac{0.5 \ln(1+1.2288 \epsilon_L)}{\left(\epsilon_L + 0.1728 \sqrt{\epsilon_L} + 0.008 \epsilon_L^{0.1504}\right)^n}\right)
!> \]
pure function calc_E_sputtered_energy_coeff(this, E) result(yield)
  class(eckstein_sputtered_energy_coeff), intent(in) :: this
  real*8, intent(in)                                 :: E !< energy in eV
  real*8 :: yield
  real*8 :: f_stopping_power
  real*8 :: eps_l !< Lindhard reduced energy
  real*8 :: w
  if (.not. allocated(this%yn)) then
    yield = 0.d0
    return
  end if

  eps_l = E/this%eps
  w = eps_l + 0.1728d0 * sqrt(eps_l) + 0.008d0 * (eps_l)**0.1504d0

  f_stopping_power = (0.5d0 * log(1.d0+1.2288d0*E/this%eps) / w**this%nu)**this%n

  yield = this%q * f_stopping_power * (E/this%E_threshold - 1.d0)**this%mu / (this%lambda/w + (E/this%E_threshold - 1.d0)**this%mu)
end function calc_E_sputtered_energy_coeff



!> Read a file containing eckstein coefficients in the following format:
!>
!> HEADER 1
!> sputter yield / sputtered energy coefficients
!> [blank]
!> HEADER 2
!> [blank]
!> [blank]
!> Angle dependent coefficient line 1
!> ...
!>
!> An example is given here
!> ```
!> Ion Target lambda q mu Eth epsilon Esb Esb/gamma
!> D W 0.3583 0.0183 1.4410 228.84 9.92326e+3 8.68 202.85
!>
!> Ion Target E0 f b c Y(E0,0) Esp theta* theta_max
!> 
!>
!> D W 250 4.2860 2.9471 0.7250 2.34e-5 1.00 93.62 44.77
!> ...
!> ```
subroutine read(this)
  use mpi_mod
  class(eckstein_coeff_set), intent(inout) :: this ! should have Z_ion and Z_target set
  character(len=2) :: ctype
  integer :: u, ierr, i
  character(len=20) :: file
  character(len=3) :: tmp
  type(eckstein_coeff), allocatable :: yn(:), temp(:) !< dummy variables to read into
  integer :: my_id

  call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)
  

  select type (c => this)
  type is (eckstein_sputter_yield)
    ctype = 'y'
  type is (eckstein_sputtered_energy_coeff)
    ctype = 'ye'
  class default
    if (my_id .eq. 0) write(*,*) 'ERROR: unknown coeff_type for eckstein coeff read'
    return
  end select


  if (allocated(this%yn)) deallocate(this%yn)

  if (this%Z_ion .lt. lbound(element_symbols,1)) then
    if (my_id .eq. 0) write(*,*) 'ERROR: unknown ions requested for eckstein coeff read', this%Z_ion, this%Z_target
    return
  end if

  write(file,"(A,A,A,A,A)") trim(ctype), "_", trim(element_symbols(this%Z_ion)), trim(element_symbols(this%Z_target)), ".dat"


  ! newunit=u breaks something in gfortran 7.2.1 combined with above internal write, so I'm doing it the oldschool way
  u = 10
  open(unit=u,file=file,status='old',access='sequential',form='formatted', action='read', iostat=ierr)
  if (ierr .ne. 0) then
    if (my_id .eq. 0) write(*,*) 'ERROR: cannot read file ', file, ' setting yields to 0'
    return
  end if

  ! Read and ignore header
  read(u,*) ! Ion Target lambda q mu Eth epsilon Esb Esb/gamma

  ! Depending on the type read the normal incidence fit
  select type (c => this)
  type is (eckstein_sputter_yield)
    read(u,*) tmp, tmp, c%lambda, c%q, c%mu, c%E_threshold, c%eps, tmp, tmp
  type is (eckstein_sputtered_energy_coeff)
    read(u,*) tmp, tmp, c%lambda, c%q, c%mu, c%E_threshold, c%nu, c%eps, tmp
    c%n = 1
    if (c%Z_ion .gt. 2) c%n = 2
  end select
  read(u,*) ! empty line

  ! Read and ignore second header
  read(u,*) ! Ion Target E0 f b c Y(E0,0) Esp theta* theta_max
  read(u,*) ! empty line
  read(u,*) ! empty line

  ! Allocate temporary storage
  allocate(yn(64))

  i = 0
  do
    i = i+1
    yn(i)%E0 = 0.d0
    read(u,*,iostat=ierr) tmp, tmp, yn(i)%E0, yn(i)%f, yn(i)%b, yn(i)%c, &
        yn(i)%Yn, yn(i)%Esp, yn(i)%theta_star, yn(i)%theta_max
    if (ierr .gt. 0) then
      write(*,*) 'ERROR: error reading file ', file
      exit
    else if (ierr .lt. 0 .or. abs(yn(i)%E0) .le. 1d-100) then ! end of file
      i = i-1
      exit
    end if
    if (i .ge. size(yn,1)) then
      allocate(temp(2*size(yn,1)))
      temp(:size(yn,1)) = yn
      call move_alloc(temp, yn) ! temp gets deallocated
    end if

    ! Ignore this line if it is a duplicate
    if (any(abs(yn(i)%E0-yn(1:i-1)%E0) .le. 1d-10)) i=i-1
  end do

  allocate(this%yn(i))
  this%yn = yn(:i)
  close(unit=u)
end subroutine read


!> Function for evaluating the angle dependency of either the sputtering yield or the sputtered energy coefficient
!> Only 1 function for the energy dependency is need, as the sputter yield and the sputter energy coefficient have the same fit formula
!> just different values for the fit parameters
pure function evaluate_eckstein_formula(this, theta) result(coeff) !< this , theta
  class(eckstein_coeff), intent(in) :: this
  real*8 ,intent(in)                :: theta !< angle in degrees
  real*8                            :: coeff
  ! Normalized angle of incidence dependency Y(E_0, theta_0)/Y(E_0,0)
  if (abs(theta) .gt. 1d-6) then
    coeff = (cos( (theta/this%theta_star * pi/2.d0)**this%c) )**(-this%f) * &
        exp(this%b * (1.d0-1.d0/cos( (theta/this%theta_star * pi/2.d0)**this%c )))
  else
    coeff = 1.d0
  end if
end function evaluate_eckstein_formula
end module mod_eckstein_y_ye
