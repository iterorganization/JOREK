!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> (M)onte-(C)arlo (C)oulomb (C)ollisions. Computes the test-particle Monte Carlo
!> coefficients and advances the related stochastic differential equations.
!> The maximum accepted normalized temperature of the background species is set to 
!> Theta=T/mc^2=0.1. This covers every  Tokamak and stellarator.
!<
module mccc_mod
  use interpN, only: interp2
  use special_mod, only : special_L0L1, special_besk0exp, special_besk1exp, special_besk2exp
  implicit none
  
  private
  real*8, parameter :: pi = 4.D0*atan(1.D0) 
  real*8, parameter :: eps0 = 8.854187817D-12 ![F/m]
  real*8, parameter :: c = 299792458.D0 ![m/s]
  real*8, parameter :: elemCharge = 1.60217662D-19 ![C]
  real*8, parameter :: hbar = 1.05457180D-34

  real*8, parameter, private :: MCCC_NANRESULT = -3
  
  type mccc_special
     real*8, allocatable, dimension(:)   :: u
     real*8, allocatable, dimension(:)   :: theta
     real*8, allocatable, dimension(:,:) :: L0
     real*8, allocatable, dimension(:,:) :: L1
  end type mccc_special

  public :: mccc_special, mccc_init, mccc_uninit, mccc_readdat, mccc_writedat, &
       mccc_push, mccc_clog, mccc_coeffs, mccc_error, mccc_testCoeff
  
contains

  !> Initializes the mccc_special struct that stores the data required for
  !> interpolating integrals L0 and L1. The interpolation grid is
  !> logarithmic and uses linear interpolation. Deallocate with mccc_uninit().
  !>
  !> input:
  !>
  !> real*8  uminxp      -- defines minimum u as umin=10^uminxp
  !> real*8  umaxxp      -- defines maximum u as umax=10^umaxxp
  !> real*8  thminxp     -- defines minimum theta as thmin=10^thminxp
  !> real*8  thmaxxp     -- defines maximum theta as thmax=10^thmaxxp
  !> integer nu          -- number of u grid points
  !> integer nth         -- number of theta grid points
  !> real*8  L0L1_eps    -- (optional) 
  !> real*8  L0L1_cutoff -- (optional) 
  !>
  !> output:
  !> 
  !> mccc_special        -- interpolation data
  !<
  type(mccc_special) function mccc_init(uminxp,umaxxp,thminxp,thmaxxp,nu,nth,&
       L0L1_eps,L0L1_cutoff)
    implicit none
    real*8, intent(in)  :: uminxp,umaxxp,thminxp,thmaxxp     
    integer, intent(in) :: nu,nth        
    real*8, intent(inout), optional :: L0L1_eps, L0L1_cutoff

    integer :: i, j
    
    allocate(mccc_init%u(nu),mccc_init%theta(nth),mccc_init%L0(nu,nth),mccc_init%L1(nu,nth))

    mccc_init%u=10**(uminxp + (/ ( (i-1)*(umaxxp-uminxp)/(nu-1), i=1,nu) /) )
    mccc_init%theta=10**(thminxp + (/ ( (i-1)*(thmaxxp-thminxp)/(nth-1), i=1,nth) /))

    do i=1,nu
       do j=1,nth
          if(present(L0L1_eps) .and. present(L0L1_cutoff)) then
             call special_L0L1(mccc_init%u(i),mccc_init%theta(j),mccc_init%L0(i,j),mccc_init%L1(i,j),&
                  L0L1_eps=L0L1_eps,L0L1_cutoff=L0L1_cutoff)
          elseif(present(L0L1_eps)) then
             call special_L0L1(mccc_init%u(i),mccc_init%theta(j),mccc_init%L0(i,j),mccc_init%L1(i,j),&
                  L0L1_eps=L0L1_eps)
          elseif(present(L0L1_cutoff)) then
             call special_L0L1(mccc_init%u(i),mccc_init%theta(j),mccc_init%L0(i,j),mccc_init%L1(i,j),&
                  L0L1_cutoff=L0L1_cutoff)
          else
             call special_L0L1(mccc_init%u(i),mccc_init%theta(j),mccc_init%L0(i,j),mccc_init%L1(i,j))
          end if
       end do
       write(*,*) "done: ",i,"/",nu
    end do
    
  end function mccc_init

  !> Deinitializes the mccc_special struct.
  !>
  !> input:
  !>
  !> mccc_special dat -- interpolation data
  !<
  integer function mccc_uninit(dat)
    implicit none
    type(mccc_special), intent(inout) :: dat
   
    deallocate(dat%u,dat%theta,dat%L0,dat%L1)

    mccc_uninit = 0
  end function mccc_uninit

  !> Writes the mccc_special struct to a file.
  !>
  !> input:
  !>
  !> mccc_special dat      -- interpolation data
  !> character(len=*) fn   -- output filename
  !>
  !> output:
  !>
  !> integer mccc_writedat -- negative if there was an error
  !<
  integer function mccc_writedat(dat,fn)
    implicit none
    type(mccc_special), intent(in)     :: dat
    character(len=*), intent(in) :: fn

    integer :: chn=989
    integer :: nu,nth

    mccc_writedat = 0

    nu = size(dat%u)
    nth = size(dat%theta)

    open(unit=chn,file=fn,action='write', err=100)
    write(chn,*) nu
    write(chn,*) nth
    write(chn,*) dat%u
    write(chn,*) dat%theta
    write(chn,*) dat%L0
    write(chn,*) dat%L1
    close(chn, err=100)

    100 mccc_writedat = -1
    
  end function mccc_writedat

  !> Reads the mccc_special struct from a file.
  !>
  !> input:
  !>
  !> character(len=*) fn -- input filename
  !>
  !> output
  !>
  !> mccc_special        -- interpolation data
  !<
  type(mccc_special) function mccc_readdat(fn)
    character(len=*), intent(in) :: fn

    integer :: chn=989
    integer :: nu,nth

    open(unit=chn,file=fn,action='read')
    read(chn,*) nu
    read(chn,*) nth
    allocate(mccc_readdat%u(nu),mccc_readdat%theta(nth),&
         mccc_readdat%L0(nu,nth),mccc_readdat%L1(nu,nth))
    read(chn,*) mccc_readdat%u
    read(chn,*) mccc_readdat%theta
    read(chn,*) mccc_readdat%L0
    read(chn,*) mccc_readdat%L1
    close(chn)
    
  end function mccc_readdat

  !> Computes Coulomb logarithm for a given test particle and plasma species.
  !> The logarithm is estimated as ln{lambda_D/min{bqm,bcl}} where lambda_D is
  !> the Debye length, and bqm and bcl are quantum mechanical and classical
  !> impact parameters, respectively. The Coulomb logarithm for different
  !> plasma species is returned.
  !>
  !> input:
  !> 
  !> mccc_special dat -- mccc_special struct that is obtained with the mccc_init() call
  !> real*8 ma        -- test particle mass [kg]
  !> real*8 qa        -- test particle charge [C]
  !> real*8 mb(:)     -- list of background species masses [kg]
  !> real*8 qb(:)     -- list of background species charges [C]
  !> real*8 nb(:)     -- list of background densities [1/m^3]
  !> real*8 thb(:)    -- list of normalized background temperatures thb=T_b/(m_b*c^2) [1] 
  !> real*8 u         -- normalized test particle momentum magnitude u=|p|/mc [1]
  !>
  !> output:
  !>
  !> real*8 clogab(:) -- list of Coulomb logarithms for species a colliding with b [1] 
  !<
  subroutine mccc_clog(ma,qa,mb,qb,nb,thb,u,clog)
    implicit none
    real*8, intent(in)  :: ma, qa, u
    real*8, intent(in)  :: mb(:), qb(:)
    real*8, intent(in)  :: thb(:), nb(:)
    real*8, intent(out) :: clog(:)

    real*8  :: debyeLength ! Debye length accounting for all plasma species
    real*8  :: mr          ! Reduced mass 
    real*8  :: bcl         ! Classical impact parameter
    real*8  :: bqm         ! Quantum mechanical impact parameter
    integer :: nspec,i     ! Helper variables
    real*8, dimension(size(clog)) :: ubar ! Mean relative velocity

    debyeLength = sqrt(eps0*c**2/sum( (nb*qb**2)/(thb*mb) ))
    nspec = size(clog)
    ubar = c*sqrt(u**2/(1+u**2)*(/(1,i=1,nspec)/)+3*thb)

    do i=1,nspec
       mr = ma*mb(i)/(ma+mb(i))
       bcl = qa*qb(i)/(4*pi*eps0*mr*ubar(i)**2)
       bqm = hbar/(2*mr*ubar(i))

       clog(i) = log(debyeLength/max(bcl,bqm))
    end do

  end subroutine mccc_clog

  !> Computes requested collision coefficients and respective derivatives
  !> for a given test particle and background species.
  !> 
  !> input:
  !> 
  !> mccc_special dat -- mccc_special struct that is obtained with the mccc_init() call
  !> real*8 ma        -- test particle mass [kg]
  !> real*8 qa        -- test particle charge [C]
  !> real*8 clogab    -- coulomb logarithm for species a colliding with b [1] 
  !> real*8 mb        -- background species mass [kg]
  !> real*8 qb        -- background species charge [C]
  !> real*8 nb        -- background density [1/m^3]
  !> real*8 thb       -- normalized background temperature thb=T_b/(m_b*c^2) [1] 
  !> real*8 u         -- normalized test particle momentum magnitude u=|p|/mc [1]
  !>
  !> output:
  !> 
  !> real*8 K         -- (optional) friction coefficient [1/s]
  !> real*8 dK        -- (optional) derivative of K with respect to u [1/s]
  !> real*8 Dpar      -- (optional) parallel momentum diffusion [1/s]
  !> real*8 dDpar     -- (optional) derivative of Dpar with respect to u [1/s]
  !> real*8 Dperp     -- (optional) perpendicular momentum diffusion [1/s]
  !> real*8 dDperp    -- (optional) derivative of Dperp with respect to u [1/s]
  !> real*8 kappa     -- (optional) guiding center friction coefficient [1/s]
  !> real*8 dkappa    -- (optional) derivative of kappa with respect to u [1/s]
  !<
  subroutine mccc_coeffs(dat,ma,qa,clogab,mb,qb,nb,thb,u,&
       K,dK,Dpar,dDpar,Dperp,dDperp,kappa,dkappa)
    implicit none
    type(mccc_special), intent(in) :: dat
    real*8, intent(in) :: ma, qa, u
    real*8, intent(in) :: mb, qb,clogab, nb, thb
    real*8, intent(out), optional :: K,Dpar,Dperp, kappa
    real*8, intent(out), optional :: dK,dDpar,dDperp, dkappa

    real*8  :: Gab,mu0,mu1,mu2,gamma,dmu0,dmu1,dmu2 ! Special functions and coefficients
    real*8  :: u2,u3,u4,gamma2,gamma3 ! Helper variables
    
    gamma=sqrt(1.D0+u**2)
    Gab=nb*(qa*qb)**2*clogab/(4*pi*eps0**2*ma**2*c**3)
    u2=u**2
    u3=u**3

    if(present(dK) .or. present(dDpar) .or. present(dDperp) .or. present(dkappa)) then
       u4=u**4
       gamma2=gamma**2
       gamma3=gamma**3
       call mccc_mufuncs(u,thb,dat,mu0,mu1,mu2,dmu0,dmu1,dmu2)
    else
       call mccc_mufuncs(u,thb,dat,mu0,mu1,mu2)
    end if
    
    if(present(K)) then
       K=-Gab*(mu0/gamma+(ma/mb)*mu1)/u2
    end if
    if(present(dK)) then
       dK=(Gab/u3)*(2*(mu0/gamma+(ma/mb)*mu1)-u*(dmu0/gamma+(ma/mb)*dmu1)+u2*mu0/gamma**3)
    end if

    if(present(Dpar)) then
       Dpar=Gab*gamma*thb*mu1/u3
    end if
    if(present(dDpar)) then
       dDpar=(Gab*thb/(gamma*u4))*(gamma2*u*dmu1-(1.D0+2*gamma2)*mu1)
    end if

    if(present(Dperp)) then
       Dperp=Gab*(u2*(mu0+gamma*thb*mu2)-thb*mu1)/(2*gamma*u3)
    end if
    if(present(dDperp)) then
       dDperp=(Gab/(2*gamma3*u4))*((4*gamma2-1.D0)*thb*mu1-u2*((2*gamma2-1.D0)*mu0+thb*&
            gamma3*mu2)+gamma2*(u3*(dmu0+thb*gamma*dmu2)-thb*u*dmu1))
    end if

    if(present(kappa)) then
       kappa=-Gab*(ma/mb)*mu1/u2
    end if
    if(present(dkappa)) then
       dkappa=-Gab*(ma/mb)*(dmu1-2*mu1/u)/u2
    end if

  end subroutine mccc_coeffs
  
  !> Computes the value for particle momentum after collisions with
  !> background species using Euler-Maruyama method with a fixed time step.
  !> 
  !> input:
  !> 
  !> mccc_special dat  -- mccc_special struct that is obtained with the mccc_init() call
  !> real*8  ma        -- test particle mass [kg]
  !> real*8  qa        -- test particle charge [C]
  !> real*8  clogab(:) -- list of coulomb logarithms for species a colliding with b [1] 
  !> real*8  mb(:)     -- list of background species masses [kg]
  !> real*8  qb(:)     -- list of background species charges [C]
  !> real*8  nb(:)     -- list of background densities [1/m^3]
  !> real*8  thb(:)    -- list of normalized background temperatures thb=T_b/(m_b*c^2) [1] 
  !> real*8  dt        -- time step length [s]
  !> real*8  rnd(3)    -- array with three elements of standard normal random numbers ~ N(0,1)
  !> real*8  uin(3)    -- normalized test particle momentum u=p/mc before collisions [1]
  !>
  !> output:
  !> 
  !> real*8  uout(3)   -- normalized test particle momentum u=p/mc after collisions [1]
  !> integer err       -- error flag, negative indicates something went wrong
  !<
  subroutine mccc_push(dat,ma,qa,clogab,mb,qb,nb,thb,dt,rnd,uin,uout,err)
    implicit none
    type(mccc_special), intent(in) :: dat 
    real*8, intent(in) :: ma, qa         
    real*8, intent(in) :: mb(:), qb(:), nb(:), thb(:), clogab(:)    
    real*8, intent(in) :: dt       
    real*8, intent(in) :: rnd(3)    
    real*8, intent(in) :: uin(3)  
    
    real*8, intent(out)  :: uout(3)
    integer, intent(out) :: err
    

    real*8 :: K,Dpar,Dperp,Kb,Dparb,Dperpb ! the fokker-planck coefficients
    real*8 :: dW(3)                        ! the change in the Wiener process during dt
    real*8 :: uhat(3)                      ! a unit vector parallel to pin
    real*8 :: u                           ! absolute value of particle momentum normalized to mc
    integer  :: i, nspecies                ! for iterating over plasma species

    err = 0

    ! Wiener process for this step
    dW=sqrt(dt)*rnd

    ! Magnitude and unit vector of p
    u=norm2(uin)
    uhat=uin/u
    
    ! Evaluate and sum Fokker-Planck coefficients
    K=0.D0
    Dpar=0.D0
    Dperp=0.D0
    nspecies = size(clogab)
    do i = 1,nspecies
       call mccc_coeffs(dat,ma,qa,clogab(i),mb(i),qb(i),nb(i),thb(i),u,K=Kb,Dpar=Dparb,Dperp=Dperpb)
       K=K+Kb
       Dpar=Dpar+Dparb
       Dperp=Dperp+Dperpb
    end do

    ! Use Euler-Maruyama method to get pout
    uout = uin+K*uhat*dt+sqrt(2*Dpar)*dot_product(uhat,dW)*uhat+sqrt(2*Dperp)*(dW-dot_product(uhat,dW)*uhat)

    if(any(isnan(uout))) err = MCCC_NANRESULT

  end subroutine mccc_push

  !> Prints error diagnosis if an error has occurred.
  !>
  !> input:
  !>
  !> integer err -- error flag, negative indicates something went wrong
  !<
  subroutine mccc_error(err)
    implicit none
    integer, intent(in) :: err

    if(err .ge. 0) then
       return
    elseif(err .eq. MCCC_NANRESULT) then
       print*, 'Error: Result is NaN.'
    else
       print*, 'Error: Unknown error.'
    end if

  end subroutine mccc_error

  !> Evaluates the mu functions
  subroutine mccc_mufuncs(u,th,data,mu0,mu1,mu2,dmu0,dmu1,dmu2)
    implicit none
    real*8, intent(in) :: u, th
    type(mccc_special), intent(in) :: data
    real*8, intent(out) :: mu0,mu1,mu2
    real*8, intent(out), optional :: dmu0,dmu1,dmu2
    real*8 :: gamma,gammasq,expBessel2,L0,L1,expgammatheta,tg,th2,u2,tgK

    gammasq=1.D0+u**2
    gamma=sqrt(gammasq)
    expBessel2=special_besk2exp(1.D0/th)
    call mccc_L0L1(u,th,data,L0,L1)
    expgammatheta=exp((1.D0-gamma)/th)
    tg=th*gamma
    th2=th**2

    mu0=(gammasq*L0-th*L1+(th-gamma)*u*expgammatheta)/expBessel2
    mu1=(gammasq*L1-th*L0+(th*gamma-1)*u*expgammatheta)/expBessel2
    mu2=(2*tg*L1+(1+2*th2)*u*expgammatheta)/(th*expBessel2)

    if(present(dmu0) .or. present(dmu1) .or. present(dmu2)) then
       u2=u**2
       tgK=tg*expBessel2
       if(present(dmu0)) then
          dmu0=(2*tg*u*L0+(gamma-2*th)*u2*expgammatheta)/tgK
       end if
       if(present(dmu1)) then
          !dmu1=(2.D0*tg*u*L1+(2.D0*th2+1.D0)*u2*expgammatheta)/tgK ! The explicit form
          dmu1=mu2*u/gamma
       end if
       if(present(dmu2)) then
          dmu2=(2*th2*u*L1+(2*th2*tg+2*th2+tg-u2)*expgammatheta)/(tgK*th)
       end if
    end if
    
  end subroutine mccc_mufuncs


  !> Evaluates the special functions L0 and L1
  !> using approximations outside the tabulated domain
  !> if theta>10**thetamaxXponent is requested, uses
  !> theta = 10**thetamaxXponent, and gives a warning
  subroutine mccc_L0L1(u,theta,data,L0,L1)
    type(mccc_special), intent(in) :: data
    real*8, intent(in) :: u,theta
    real*8, intent(out) :: L0,L1
    real*8 :: th
    integer :: nu, nth

    nu = size(data%u)
    nth = size(data%theta)

    if ( theta.gt.data%theta(nth) ) then
       th=data%theta(nth) 
       print*,'Warning, temperature exceeds tabulated'
    else
       th=theta
    end if
    
    if(u.gt.data%u(nu)) then
       L0=special_besk0exp(1.D0/th)
       L1=special_besk1exp(1.D0/th)
    else
       if( (u.gt.data%u(1)) .and. (th.gt.data%theta(1)) ) then
          L0 = interp2(data%u,data%theta,data%L0,u,th,nu,nth) 
          L1 = interp2(data%u,data%theta,data%L1,u,th,nu,nth) 
       else
          L0 = sqrt(pi*theta/2)*erf(u/sqrt(2*theta))
          L1=L0
       end if
    end if
    
  end subroutine mccc_L0L1

  !> This test program evaluates Fokker-Planck coefficients and respective derivatives
  !> (analytically) for some theta and a range of u values.
  !> theta is background normalized temperature T_b/(m_b*c^2) and u is normalized
  !> test particle momentum p/mc. The results are then written in mccc_coefficients.out
  !<
  subroutine mccc_testCoeff
    implicit none

    real*8, parameter :: uminxp = -3.D0
    real*8, parameter :: umaxxp = 2.D0
    real*8, parameter :: thminxp = -3.D0
    real*8, parameter :: thmaxxp = -1.D0 
    integer, parameter :: nu = 401
    integer, parameter :: nth = 201

    integer, parameter :: nspec=1
    type(mccc_special) :: dat    
    real*8 :: clogab(nspec) ! [1]
    real*8 :: mb(nspec) ! [kg]
    real*8 :: qb(nspec) ! [C]
    real*8 :: nb(nspec) ! [1/m^3]
    real*8 :: thb(nspec) ! [1]
    real*8 :: m ! [kg]
    real*8 :: q ! [C}
    real*8 :: dx,xmax,xmin! [p/mc]

    ! for input output
    logical :: storage_file_on_disk
    integer, parameter :: chn=989
    character(20), parameter :: storage_file='mccc.data'
    
    
    real*8 :: K,dK,Dpar,dDpar,Dperp,dDperp,kappa,dkappa,Kb,dKb,Dparb,dDparb,Dperpb,dDperpb,kappab,dkappab

    integer Niter,i,j,err

    print*,''
    print*,'This test program evaluates Fokker-Planck coefficients and respective derivatives'
    print*,'for a range of u values with a given theta. The results are written in test_coefficients.out'
    print*,''

    ! Define the background distribution to be electrons
    ! at temperature T/mc^2 = 0.1
    clogab=18.D0
    mb = 9.10938356D-31
    qb = -elemCharge
    nb = 1.D20
    thb = 1.D-1

    m = 9.10938356D-31
    q = -elemCharge
    
    xmax = 1.D1
    xmin = 1.D-3
    Niter = 1.D5
    dx = (xmax-xmin)/Niter

    ! Initialize the look-up tables
    print*,''
    write(*,*) 'Initializing look-up tables...'
    inquire(file=storage_file,exist=storage_file_on_disk)
    if(storage_file_on_disk) then
       dat = mccc_readdat(storage_file)
    else
       dat = mccc_init(uminxp,umaxxp,thminxp,thmaxxp,nu,nth)
       err = mccc_writedat(dat,storage_file)
    end if
    print*,'Done!'

    open(999,file="test_coefficients.out")
    
    
    !write(999,'(10A23)') 'u','theta','K','dK','Dpar','dDpar','Dperp','dDper','kappa','dkappa'
    
    do j=1,niter
       K=0.D0
       Dpar=0.D0
       Dperp=0.D0
       dK=0.D0
       dDpar=0.D0
       dDperp=0.D0
       kappa=0.D0
       dkappa=0.D0
       do i = 1,nspec
          call mccc_coeffs(dat,m,q,clogab(i),mb(i),qb(i),nb(i),thb(i),xmin+j*dx,&
               K=Kb,Dpar=Dparb,Dperp=Dperpb,dK=dKb,dDpar=dDparb,dDperp=dDperpb,&
               kappa=kappab,dkappa=dkappab)
          K=K+Kb
          Dpar=Dpar+Dparb
          Dperp=Dperp+Dperpb
          dK=dK+dKb
          dDpar=dDpar+dDparb
          dDperp=dDperp+dDperpb
          kappa=kappa+kappab
          dkappa=dkappa+dkappab
       end do

       write(999,'(10E23.12)') xmin+j*dx,sum(thb),K,dK,Dpar,dDpar,Dperp,dDperp,kappa,dkappa

    end do
    close(999)
    err = mccc_uninit(dat)
    
  end subroutine mccc_testCoeff

end module mccc_mod
