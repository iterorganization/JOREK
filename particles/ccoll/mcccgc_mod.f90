!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> MCCC operator for guiding center coordinates.
!> These coordinates are location X, momentum p (scalar), and
!> pitch = p_para/p. The collisions can be evaluated separately
!> for each coordinate. Both fixed and adaptive scheme are supported,
!> and Milstein method is used for both.
!<
module mcccgc_mod
  use mccc_mod
  implicit none

  real*8, parameter, private :: pi = 4.D0*atan(1.D0) 
  real*8, parameter, private :: eps0 = 8.854187817D-12 ![F/m]
  real*8, parameter, private :: c = 299792458.D0 ![m/s]
  real*8, parameter, private :: elemCharge = 1.60217662D-19 ![C]

  integer, parameter, private :: MCCCGC_NANRESULT        = -3
  integer, parameter, private :: MCCCGC_EXTREMELYSMALLDT = -4
  integer, parameter, private :: MCCCGC_NOTPHYSICAL     = -5

  real*8, parameter, private :: MCCCGC_EXTREMELYSMALLDTVAL = 1.D-12

  type mcccgc_coefstruct 
     real*8 :: kappa  = 0.D0
     real*8 :: dkappa = 0.D0
     real*8 :: Dpar   = 0.D0
     real*8 :: dDpar  = 0.D0
     real*8 :: Dperp  = 0.D0
     real*8 :: dDperp = 0.D0
     real*8 :: gx     = 0.D0
     real*8 :: nu     = 0.D0
  end type mcccgc_coefstruct

  public :: mcccgc_push, mcccgc_evalCoefs

contains

  !> Computes the collision coefficients and stores them in a struct that is used
  !> by mcccgc_push(). The coefficients are calculated separately like this so that
  !> they are not unnecessarily re-calculated if the time step is rejected.
  !> 
  !> input:
  !> 
  !> mccc_special dat -- struct that is obtained with the mccc_init() call
  !> real*8 ma        -- test particle mass [kg]
  !> real*8 qa        -- test particle charge [C]
  !> real*8 clogab(:) -- list of coulomb logarithms for species a colliding with b [1] 
  !> real*8 mb(:)     -- list of background species masses [kg]
  !> real*8 qb(:)     -- list of background species charges [C]
  !> real*8 nb(:)     -- list of background densities [1/m^3]
  !> real*8 thb(:)    -- list of normalized background temperatures thb=T_b/(m_b*c^2) [1] 
  !> real*8 u         -- normalized guiding center momentum magnitude u=p/mc [1]
  !> real*8 xi        -- guiding center pitch [1]
  !>
  !> output:
  !> 
  !> mcccgc_coefstruct coefs -- struct containing the coefficients needed for evaluating the collisions
  !<
  subroutine mcccgc_evalCoefs(dat,ma,qa,clogab,mb,qb,nb,thb,u,xi,coefs)
    type(mccc_special), intent(in) :: dat
    real*8, intent(in) :: ma,qa,u,xi
    real*8, intent(in) :: clogab(:),mb(:),qb(:),nb(:),thb(:) 

    type(mcccgc_coefstruct), intent(out) :: coefs

    real*8 :: kappa,dkappa,Dpar,dDpar,Dperp,dDperp
    real*8 :: kappab,dkappab,Dparb,dDparb,Dperpb,dDperpb

    integer :: nspecies,i

    nspecies = size(clogab)
    kappa=0.D0
    dkappa=0.D0
    Dpar=0.D0
    dDpar=0.D0
    Dperp=0.D0
    dDperp=0.D0
    do i = 1,nspecies
       call mccc_coeffs(dat,ma,qa,clogab(i),mb(i),qb(i),nb(i),thb(i),u,kappa=kappab,Dpar=Dparb,Dperp=Dperpb,&
         dkappa=dkappab,dDpar=dDparb,dDperp=dDperpb)
       
       kappa=kappa+kappab
       dkappa=dkappa+dkappab
       Dpar=Dpar+Dparb
       dDpar=dDpar+dDparb
       Dperp=Dperp+Dperpb
       dDperp=dDperp+dDperpb
       
    end do
    coefs%kappa  = kappa
    coefs%dkappa = dkappa
    coefs%Dpar   = Dpar
    coefs%dDpar  = dDpar
    coefs%Dperp  = Dperp
    coefs%dDperp = dDperp

    coefs%gx = (ma*c)**2*( (Dpar-Dperp)*(1-xi**2)/2 + Dperp )/qa**2
    coefs%nu = 2*Dperp/u**2
  end subroutine mcccgc_evalCoefs
  
  !> Computes the value for guiding center phase space position after collisions with
  !> background species using Milstein method with a fixed or an adaptive time step.
  !> Returns a suggestion for the next time step. A negative suggestion indicates 
  !> the current time step was rejected, and absolute value should be used for the 
  !> next try. If no tolerance is given, Milstein method with fixed time step will
  !> be used instead. Wiener processes are generated automatically and only cleaning
  !> is required after each accepted time step.
  !>
  !> input:
  !> 
  !> mcccgc_coefstruct coefs -- struct containing the coefficients from mccca_evalCoefs() call
  !> real*8 wienarr(7,:)     -- array that stores the Wiener processes (Ndim=5)
  !> real*8 B(3)             -- magnetic field at guiding center position before collisions [T]
  !> real*8 cutoff           -- value (>= 0) used to mirror u as "u = 2*cutoff - u" if u < cutoff otherwise [1]
  !> real*8 uin              -- normalized test particle momentum u=p/mc before collisions [1]
  !> real*8 xiin             -- guiding center pitch before collisions [1]
  !> real*8 dtin             -- candidate time step length [s]
  !> real*8 tol(2)           -- (optional) two-element array containing the relative tolerance for the
  !>                            guiding center momentum and absolute tolerance for the pitch. Fixed time 
  !>                            step is used if not given
  !>
  !> output:
  !> 
  !> real*8 uout  -- normalized guiding center momentum u=p/mc after collisions [1]
  !> real*8 xiout -- guiding center pitch after collisions [1]
  !> real*8 dx(3) -- change in guiding center position due to collisions [m]
  !> real*8 dtout -- suggestion for the next time step. dtout = dtin in case fixed step is used [s]
  !> integer err  -- error flag, negative value indicates error
  !<
  subroutine mcccgc_push(coefs,cutoff,uin,uout,xiin,xiout,dtin,rnd)
    type(mcccgc_coefstruct), intent(in) :: coefs
    
    real*8, intent(in) :: dtin ! time step length [s]
    real*8, intent(in) :: uin ! p/mc test particle momentum vector (px,py,pz) before collision, normalized to mc 
    real*8, intent(in) :: xiin !
    real*8, intent(in) :: cutoff
    real*8, intent(in) :: rnd
    
    real*8, intent(out) :: uout 
    real*8, intent(out) :: xiout

    real*8 :: time
    real*8 :: dW(5),dW2(5) ! the change in the Wiener process during dt
    real*8 :: bhat(3)
    real*8 :: F, gu, gxi
    real*8 :: kappa_k, kappa_d(2), erru, gx
    real*8 :: dWopt(2),dti,alpha
    integer :: i, windex,tindex, ki,kmax
    logical rejected

    F = coefs%kappa+coefs%dDpar+2*coefs%Dpar/uin
    gu = sqrt(2*coefs%Dpar)
    gxi = sqrt((1-xiin**2)*coefs%nu)
    
    uout = uin + F*dtin !+ gu*sqrt(dtin)*rnd(1)
    xiout = xiin - xiin*coefs%nu*dtin + gxi*sqrt(dtin)*rnd

    ! Enforce boundary conditions
    if(uout .lt. cutoff) then
       uout = 2*cutoff-uout
    end if

    if(abs(xiout) .gt. 1.D0) then
       xiout = sign(2.D0-abs(xiout),xiout)
    end if


  end subroutine mcccgc_push

  !> Prints error diagnosis if an error has occurred.
  !>
  !> input:
  !>
  !> integer err -- error flag, negative indicates something went wrong
  !<
  subroutine mcccgc_error(err)
    implicit none
    integer, intent(in) :: err

    if(err .ge. 0) then
       return
    elseif(err .eq. MCCCGC_NANRESULT) then
       print*, 'Error: Result is NaN.'
    elseif(err .eq. MCCCGC_EXTREMELYSMALLDT) then
       print*, 'Error: Extremely small timestep. Below ', MCCCGC_EXTREMELYSMALLDTVAL, ' s'
    elseif(err .eq. MCCCGC_NOTPHYSICAL) then
       print*, 'Error: Result is not physical. Either |xi| > 0 or u < 0.'
    end if

  end subroutine mcccgc_error

end module mcccgc_mod
