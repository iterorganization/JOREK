!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> Coulomb collisions for (relativistic) test particles.
!> See Sarkimaki et al, "Adaptive time-stepping Monte Carlo integration of Coulomb collisions",
!> Comp. Phys. Comm.
!> For the partial screening operator, see PhD thesis by Linnea Hesslow.
!<
module mod_ccoll_relativistic
  use data_structure
  use constants
  use hdf5_io_module
  use mod_particle_types
  use mod_fields
  use mod_bessel, only : bessel_k2exp, bessel_k1exp, bessel_k0exp
  use mod_simpson, only : simpson_adaptive, func_real8_1D
  use mod_interp_methods, only: interp_bilinear
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use, intrinsic :: ieee_arithmetic
  implicit none

  real*8, parameter :: DEFAULT_L0L1_eps    = 1.D-8 !< default tolerance in eval_L0L1
  real*8, parameter :: DEFAULT_L0L1_cutoff = 1.D-7 !< default cutoff in evalL0L1
  real*8, parameter :: DEFAULT_UCUTOFF     = 1.D-4 !< Minimum p/mc value test particle can have
  
  ! Struct for storing tabulated values of special functions L0 and L1 as well as the ion data
  type ccoll_data
     real*8, allocatable, dimension(:)    :: u     !< p/mc abscissa for tabulated data
     real*8, allocatable, dimension(:)    :: theta !< T/mc^2 abscissa for tabulated data
     real*8, allocatable, dimension(:,:)  :: L0    !< Tabulated values for L0
     real*8, allocatable, dimension(:,:)  :: L1    !< Tabulated values for L1
     real*8, allocatable, dimension(:)    :: mi    !< Ion species mass
     real*8, allocatable, dimension(:)    :: Ii    !< Mean excitation energy for ions
     integer, allocatable, dimension(:)   :: ai    !< Ion normalized effective length sclae
     integer*1, allocatable, dimension(:) :: Zi    !< Ion atomic number
     integer*1, allocatable, dimension(:) :: Z0    !< Ion charge state
     real*8 :: m_i_over_m_imp                      !< Main ion mass / impurity species mass
  end type ccoll_data

  private

  public :: ccoll_data, ccoll_compute_L0L1table, ccoll_write_L0L1table, ccoll_read_L0L1table, &
       ccoll_init, ccoll_deallocate, ccoll_kinetic_relativistic_push, ccoll_gc_relativistic_push, &
       ccoll_kinetic_relativistic_explicitpush, ccoll_gc_relativistic_explicitpush, &
       ccoll_explicitpush_partialscreening, ccoll_gc_relativistic_push_partialscreening, &
       ccoll_kinetic_relativistic_push_partialscreening

contains

  !> Initializes data for collision evaluation
  !> Parameters and tabulated values are stored in the returned struct.
  subroutine ccoll_init(fn, dat)
    implicit none
    character(len=*), intent(in)  :: fn         !< Filename where L0L1 are tabulated
    type(ccoll_data), intent(inout) :: dat      !< Contains parameters and tabulated values

    dat = ccoll_read_L0L1table(fn)
    call ccoll_init_ions(dat)
  end subroutine ccoll_init

  !> Evaluates the special functions L0 and L1.
  !> The cutoff parameter (< 1) divides the integral in two
  !> parts at u_cutoff = sqrt((1-theta * ln(cutoff))**2 -1).
  !> This helps to ensure that adaptive integration does not
  !> fail when theta is very small.
  subroutine eval_L0L1(u,theta,L0,L1,eps,cutoff)
    implicit none
    real*8, intent(in) :: u                !< p/mc value where functions are evaluated
    real*8, intent(in) :: theta            !< T/mc^2 value where functions are evaluated
    real*8, intent(in), optional :: eps    !< error tolerance for the adaptive simpsons rule
    real*8, intent(in), optional :: cutoff !< cutoff value

    real*8, intent(out) :: L0 !< evaluated L0 integral (Eq. 17 in referenced paper)
    real*8, intent(out) :: L1 !< evaluated L1 integral (Eq. 18 in referenced paper)

    real*8 :: def_cutoff, u_cutoff
    real*8 :: def_tol, tol

    def_tol    = DEFAULT_L0L1_eps
    def_cutoff = DEFAULT_L0L1_cutoff
    if(present(eps))    def_tol    = eps
    if(present(cutoff)) def_cutoff = cutoff

    u_cutoff = sqrt((1.D0 - theta * log(def_cutoff))**2 - 1.D0)

    if(u_cutoff > u) then
       tol = 0.5D0 * def_tol * u
       L0 = simpson_adaptive(L0_integrand,0.D0,u,tol,10)
       L1 = simpson_adaptive(L1_integrand,0.D0,u,tol,10)
    else
       tol = 0.5D0 * def_tol * u_cutoff
       L0 = simpson_adaptive(L0_integrand,0.D0,u_cutoff,tol,10) + simpson_adaptive(L0_integrand,u_cutoff,u,tol,20)
       L1 = simpson_adaptive(L1_integrand,0.D0,u_cutoff,tol,10) + simpson_adaptive(L1_integrand,u_cutoff,u,tol,20)
    end if

  contains

    function L0_integrand(u) result (val)
      real*8, intent(in) :: u
      real*8 :: val,gamma
      gamma = sqrt( 1.D0 + u**2 )
      val   = exp( ( 1.D0 - gamma ) / theta ) / gamma
    end function L0_integrand

    function L1_integrand(u) result (val)
      real*8, intent(in) :: u
      real*8 :: val,gamma
      gamma = sqrt( 1.D0 + u**2 )
      val   = exp( (1.D0 - gamma ) / theta )
    end function L1_integrand

  end subroutine eval_L0L1

  !> Allocates and initializes tables containing computed  L0 and L1 values for interpolation
  function ccoll_compute_L0L1table(uminxp,umaxxp,thminxp,thmaxxp,nu,nth,&
       eps,cutoff) result(tabulatedL0L1) 
    implicit none

    real*8, intent(in)  :: uminxp  !> minimum u as umin=10^uminxp 
    real*8, intent(in)  :: umaxxp  !> maximum u as umax=10^umaxxp
    real*8, intent(in)  :: thminxp !> minimum theta as thmin=10^thminxp
    real*8, intent(in)  :: thmaxxp !> maximum theta as thmax=10^thmaxxp
    integer, intent(in) :: nu      !> number of u grid points
    integer, intent(in) :: nth     !> number of theta grid points
    real*8, intent(in), optional :: eps       !> error tolerance for evaluating L0 and L1
    real*8, intent(in), optional :: cutoff    !> cutoff value for evaluating L0 and L1
    real*8 :: eps_in    = DEFAULT_L0L1_eps    !> used to define default value for eps
    real*8 :: cutoff_in = DEFAULT_L0L1_cutoff !> used to define default value for cutoff

    type(ccoll_data) :: tabulatedL0L1

    integer :: i, j
   
    ! Set optional parameters 
    if(present(eps))    eps_in    = eps
    if(present(cutoff)) cutoff_in = cutoff

 
    allocate(tabulatedL0L1%u(nu), &
         tabulatedL0L1%theta(nth),&
         tabulatedL0L1%L0(nu,nth),&
         tabulatedL0L1%L1(nu,nth))

    ! Set abscissae
    tabulatedL0L1%u     = 10**( uminxp  + (/ ( ( i - 1 ) * ( umaxxp  - uminxp  ) / ( nu  - 1 ), i=1,nu)  /) )
    tabulatedL0L1%theta = 10**( thminxp + (/ ( ( i - 1 ) * ( thmaxxp - thminxp ) / ( nth - 1 ), i=1,nth) /) )

    ! Evaluate and store values to the table

    !$omp parallel default(shared) &
    !$omp private(i,j)
    !$omp single
    do i=1,nu
       do j=1,nth
         !$omp task
         call eval_L0L1(tabulatedL0L1%u(i),tabulatedL0L1%theta(j),&
              tabulatedL0L1%L0(i,j),tabulatedL0L1%L1(i,j),&
              eps=eps_in,cutoff=cutoff_in)
         !$omp end task
       end do 
    end do
    !$omp end single
    !$omp taskwait
    !$omp end parallel  
 
  end function ccoll_compute_L0L1table


  !> Deinitializes tabulated L0L1 and ion values struct.
  subroutine ccoll_deallocate(dat)
    implicit none
    type(ccoll_data), intent(inout) :: dat ! data to be deinitialized
  
    if (allocated(dat%u))     deallocate(dat%u) 
    if (allocated(dat%theta)) deallocate(dat%theta) 
    if (allocated(dat%L0))    deallocate(dat%L0) 
    if (allocated(dat%L1))    deallocate(dat%L1) 
    if (allocated(dat%mi))    deallocate(dat%mi) 
    if (allocated(dat%Z0))    deallocate(dat%Z0) 
    if (allocated(dat%Zi))    deallocate(dat%Zi) 
    if (allocated(dat%Ii))    deallocate(dat%ai) 

  end subroutine ccoll_deallocate


  !> Writes tabulated L0L1 values to a file.
  subroutine ccoll_write_L0L1table(dat,fn)
    implicit none
    class(ccoll_data), intent(in) :: dat ! data to be written
    character(len=*), intent(in)  :: fn  ! output filename

    integer(HID_T) :: file_id
    integer :: nu, nth, ierr

    nu  = size(dat%u)
    nth = size(dat%theta)

    call HDF5_create(fn, file_id, ierr)
    if( ierr .ne. 0 ) then
       write(*,*) "Could not store L0 and L1 integrals on file."
       return
    end if

    call HDF5_integer_saving(file_id,  nu,  "nu")
    call HDF5_integer_saving(file_id, nth, "nth")
    call HDF5_array1D_saving(file_id,     dat%u,  nu,       "u")
    call HDF5_array1D_saving(file_id, dat%theta, nth,      "th")
    call HDF5_array2D_saving(file_id,    dat%L0,  nu, nth, "L0")
    call HDF5_array2D_saving(file_id,    dat%L1,  nu, nth, "L1")

    call HDF5_close(file_id)
    
  end subroutine ccoll_write_L0L1table

  !> Reads tabulated L0L1 values from a file.
  type(ccoll_data) function ccoll_read_L0L1table(fn)
    character(len=*), intent(in) :: fn !< input filename

    integer(HID_T) :: file_id
    integer :: nu, nth, ierr

    call HDF5_open(fn, file_id, ierr)
    if( ierr .ne. 0 ) then
       write(*,*) "Could not read L0 and L1 integrals from the file."
    end if
    call HDF5_integer_reading(file_id,  nu,  "nu")
    call HDF5_integer_reading(file_id, nth, "nth")
    allocate(ccoll_read_L0L1table%u(nu), ccoll_read_L0L1table%theta(nth),&
         ccoll_read_L0L1table%L0(nu,nth), ccoll_read_L0L1table%L1(nu,nth))

    call HDF5_array1D_reading(file_id, ccoll_read_L0L1table%u,      "u")
    call HDF5_array1D_reading(file_id, ccoll_read_L0L1table%theta, "th")
    call HDF5_array2D_reading(file_id, ccoll_read_L0L1table%L0,    "L0")
    call HDF5_array2D_reading(file_id, ccoll_read_L0L1table%L1,    "L1")
    call HDF5_close(file_id)
    
  end function ccoll_read_L0L1table


  !> Initialize ion data
  !> This routine allocates and initializes data needed to include ions (including impurities) to the collision operator.
  !> Call this once before calling fields%calc_njTj.
  !> All initialized arrays have size Nion and the order is (/main ion, imp_0, imp_+1, imp_+2, .../)
  subroutine ccoll_init_ions(dat)
    use phys_module, only: central_mass, imp_type
    implicit none
    type(ccoll_data), intent(inout) :: dat

    ! Mean excitation energy [eV] for all charge states from neutral to fully ionized (for which we used a dummy value as it is not used in the computation)
    real*8, parameter  :: Iconst_Ar(19) = (/188.5, 219.4, 253.8, 293.4, 339.1, 394.5, 463.4, 568.0, 728.0, 795.9, 879.8, 989.9, 1138.1, 1369.5, 1791.2, 2497.0, 4677.2, 4838.2, 1.0/)
    real*8, parameter  :: Iconst_Ne(11) = (/137.2, 165.2, 196.9, 235.2, 282.8, 352.6, 475.0, 696.8,  1409.2, 1498.4, 1.0/)

    ! Normalized effective length scale for all charge states from neutral to fully ionized (for which we used a dummy value as it is not used in the computation)
    integer, parameter :: aconst_Ar(19) = (/96, 90, 84, 78, 72, 65, 59, 53, 47, 44, 41, 38, 35, 32, 27, 21, 13, 13, 1/)
    integer, parameter :: aconst_Ne(11) = (/111, 100, 90, 80, 71, 62, 52, 40, 24, 23, 1/)

    integer*1 :: atomnum_imp, i, nions

    if(with_impurities) then
       if( trim(imp_type(1)) .eq. 'Ne') then
          atomnum_imp = 10
          dat%m_i_over_m_imp = central_mass * 0.05

          allocate( dat%Ii(atomnum_imp + 2), dat%ai(atomnum_imp + 2) )
          dat%Ii(2:atomnum_imp) = Iconst_Ne
          dat%ai(2:atomnum_imp) = aconst_Ne
       elseif( trim(imp_type(1)) .eq. 'Ar') then
          atomnum_imp = 18
          dat%m_i_over_m_imp = central_mass * 0.025

          allocate( dat%Ii(atomnum_imp + 2), dat%ai(atomnum_imp + 2) )
          dat%Ii(2:atomnum_imp) = Iconst_Ar
          dat%ai(2:atomnum_imp) = aconst_Ar
       else
          ! Unknown impurity
          write(*,*) "Unknown impurity species"
          return
       end if

       nions = size(dat%Ii)
       allocate( dat%Z0(nions), dat%Zi(nions), dat%mi(nions) )
       dat%Zi(2:nions) = atomnum_imp
       dat%Z0(2:nions) = (/ (i, i=0,atomnum_imp, 1) /)
       dat%mi(2:nions) = ATOMIC_MASS_UNIT * central_mass / dat%m_i_over_m_imp

    else
       allocate( dat%Z0(1), dat%Zi(1), dat%mi(1) )
       allocate( dat%Ii(1), dat%ai(1) )
       dat%m_i_over_m_imp = 1

    end if

    dat%mi(1) = central_mass * ATOMIC_MASS_UNIT
    dat%Zi(1) = 1
    dat%Z0(1) = 1
    dat%Ii(1) = 1.0
    dat%ai(1) = 1

    dat%Ii = dat%Ii*EL_CHG ! Convert to Joules
  
  end subroutine ccoll_init_ions


  !> Computes Coulomb logarithm for a given test particle and plasma species.
  !> The logarithm is estimated as ln{lambda_D/min{bqm,bcl}} where lambda_D is
  !> the Debye length, and bqm and bcl are quantum mechanical and classical
  !> impact parameters, respectively. The Coulomb logarithm for different
  !> plasma species is returned. 
  subroutine ccoll_clog(ma,qa,mi,qi,ne,the,ni,thi,u,cloge,clogi)
    implicit none
    real*8, intent(in)  :: ma       !< test particle mass [kg]
    real*8, intent(in)  :: qa       !< test particle charge [C]
    real*8, intent(in)  :: mi(:)    !< list of background species masses [kg]
    real*8, intent(in)  :: qi(:)    !< list of background species charges [C]
    real*8, intent(in)  :: ne       !< electron density [1/m^3]
    real*8, intent(in)  :: the      !< normalized electron temperature [T_b/(m_b*c^2)]
    real*8, intent(in)  :: ni(:)    !< list of background densities [1/m^3]
    real*8, intent(in)  :: thi(:)   !< list of normalized ion temperatures [T_b/(m_b*c^2)]
    real*8, intent(in)  :: u        !< normalized test particle momentum [p/mc]
    real*8, intent(out) :: cloge    !< Coulomb logarithm for electrons
    real*8, intent(out) :: clogi(:) !< Coulomb logarithm for each ion species

    real*8  :: debyeLength ! Debye length
    real*8  :: mr          ! Reduced mass 
    real*8  :: bcl         ! Classical impact parameter
    real*8  :: bqm         ! Quantum mechanical impact parameter
    integer :: i           ! Helper variables
    real*8  :: ubar        ! Mean relative velocity

    debyeLength = sqrt( EPS_ZERO * SPEED_OF_LIGHT**2 / ( ( ne * EL_CHG**2 ) / ( the * MASS_ELECTRON ) ) )
    ubar  = SPEED_OF_LIGHT * sqrt( u**2 / ( 1 + u**2 )  + 3.d0 * the )
    mr    = ma * MASS_ELECTRON / ( ma + MASS_ELECTRON )
    bcl   = qa * EL_CHG / ( 4.d0 * PI * EPS_ZERO * mr * ubar**2 )
    bqm   = HBAR / ( 2.d0 * mr * ubar )
    cloge = log(debyeLength/max(bcl,bqm))

    do i=1,size(mi)
       ubar  = SPEED_OF_LIGHT * sqrt( u**2 / ( 1 + u**2 )  + 3.d0 * thi(i) )
       mr  = ma * mi(i) / ( ma + mi(i) )
       bcl = qa * qi(i) / ( 4.d0 * PI * EPS_ZERO * mr * ubar**2 )
       bqm = HBAR / ( 2.d0 * mr * ubar )

       clogi(i) = log(debyeLength/max(bcl,bqm))
    end do

  end subroutine ccoll_clog


  !> Computes requested collision coefficients and respective derivatives
  !> for a given test particle and background species.
  subroutine ccoll_coeffs(dat,ma,qa,clog,mb,qb,nb,thb,u,&
       K,dK,Dpar,dDpar,Dperp,dDperp,kappa,dkappa)
    implicit none

    class(ccoll_data), intent(in) :: dat !< tabulated L0L1 values
    real*8, intent(in)  :: ma   !< test particle mass [kg]
    real*8, intent(in)  :: qa   !< test particle charge [C]
    real*8, intent(in)  :: clog !< Coulomb logarithm
    real*8, intent(in)  :: mb   !< background species mass [kg]
    real*8, intent(in)  :: qb   !< background species charge [C]
    real*8, intent(in)  :: nb   !< background density [1/m^3]
    real*8, intent(in)  :: thb  !< normalized background temperature [T_b/(m_b*c^2)]
    real*8, intent(in)  :: u    !< normalized test particle momentum [p/mc]

    real*8, intent(out), optional :: K      !< friction coefficient [1/s]
    real*8, intent(out), optional :: Dpar   !< parallel momentum diffusion [1/s]
    real*8, intent(out), optional :: Dperp  !< perpendicular momentum diffusion [1/s]
    real*8, intent(out), optional :: kappa  !< guiding center friction coefficient [1/s]
    real*8, intent(out), optional :: dK     !< derivative of K with respect to u [1/s]
    real*8, intent(out), optional :: dDpar  !< derivative of Dpar with respect to u [1/s]
    real*8, intent(out), optional :: dDperp !< derivative of Dperp with respect to u [1/s]
    real*8, intent(out), optional :: dkappa !< derivative of kappa with respect to u [1/s]

    real*8  :: Gab,mu0,mu1,mu2,gamma,dmu0,dmu1,dmu2 ! Special functions and coefficients
    real*8  :: u2,u3,u4,gamma2,gamma3 ! Helper variables
    
    gamma = sqrt(1.D0+u**2)
    Gab = nb * ( qa * qb )**2 * clog / ( 4.d0 * pi * EPS_ZERO**2 * ma**2 * SPEED_OF_LIGHT**3 )
    u2  = u**2
    u3  = u**3

    if(present(dK) .or. present(dDpar) .or. present(dDperp) .or. present(dkappa)) then
       u4 = u**4
       gamma2 = gamma**2
       gamma3 = gamma**3
       call ccoll_mufuncs(dat,u,thb,mu0,mu1,mu2,dmu0,dmu1,dmu2)
    else
       call ccoll_mufuncs(dat,u,thb,mu0,mu1,mu2)
    end if
    
    if(present(K)) then
       K = -Gab * ( mu0 / gamma + ( ma / mb) * mu1 ) / u2
    end if

    if(present(Dpar)) then
       Dpar = Gab * gamma * thb * mu1 / u3
    end if
    if(present(dDpar)) then
       dDpar = ( Gab * thb / ( gamma * u4 ) ) * ( gamma2 * u * dmu1 - ( 1.D0 + 2.d0 * gamma2 ) * mu1 )
    end if

    if(present(Dperp)) then
       Dperp = Gab * ( u2 * ( mu0 + gamma * thb * mu2 ) - thb * mu1 ) / ( 2.d0 * gamma * u3 )
    end if
    if(present(dDperp)) then
       dDperp = ( Gab / ( 2.d0 * gamma3 * u4 ) ) * ( ( 4.d0 * gamma2 - 1.D0 ) * thb * mu1 &
            - u2 * ( ( 2.d0 * gamma2 - 1.D0 ) * mu0 + thb * gamma3 * mu2 ) &
            + gamma2 * ( u3 * ( dmu0 + thb * gamma * dmu2 ) - thb * u * dmu1 ) )
    end if

    if(present(kappa)) then
       kappa = -Gab * ( ma / mb ) * mu1 / u2
    end if
    if(present(dkappa)) then
       dkappa = -Gab * ( ma / mb ) * ( dmu1 - 2.d0 * mu1 / u ) / u2
    end if

  end subroutine ccoll_coeffs


  !> Updates particle momentum after collisions
  !> Pushing is done by calling the explicit push function. This function is just a wrapper
  !> that additionally evaluates the plasma quantities and takes care of the coordinate transformation
  !> in momentum space thus simplifying the process of including collisions in simulations.
  subroutine ccoll_kinetic_relativistic_push(dat, prt, fields, mass, time, dt)
    implicit none
    class(ccoll_data), intent(in) :: dat !< Collision data
    class(particle_kinetic_relativistic), intent(inout) :: prt
    class(fields_base), intent(in) :: fields
    real*8,intent(in) :: mass, time, dt !< Mass in AMU and time in seconds

    real*8 :: E(3), B(3), psi, U, ne, rnd(3), pout(3), Te, Ti, the
    real*8, allocatable :: ni(:), thi(:)

    allocate(ni(size(dat%mi)), thi(size(dat%mi)))
    call fields%calc_NjTj(time, prt%i_elm, prt%st, prt%x(3), dat%m_i_over_m_imp, ne, Te, ni, Ti)
    the = Te * K_BOLTZ / ( MASS_ELECTRON * SPEED_OF_LIGHT**2 )
    thi = Ti * K_BOLTZ / ( dat%mi * SPEED_OF_LIGHT**2 )
    
    ! This should be ~N(0,1) but this approximation works well enough
    call random_number(rnd)
    rnd = floor(2.d0*rnd)
    rnd = -1.d0 + 2.d0 * rnd

    call ccoll_kinetic_relativistic_explicitpush(dat, mass * ATOMIC_MASS_UNIT, prt%q, &
         ne, the, ni, thi, dt, rnd, prt%p / (mass * SPEED_OF_LIGHT), pout)
    prt%p = pout * (mass * SPEED_OF_LIGHT)
    deallocate(ni, thi)
    
  end subroutine ccoll_kinetic_relativistic_push
  
  !> Computes the value for particle momentum after collisions with
  !> background species using Euler-Maruyama method with a fixed time step.
  subroutine ccoll_kinetic_relativistic_explicitpush(dat,ma,qa,ne,the,ni,thi,dt,rnd,uin,uout)
    implicit none
    class(ccoll_data), intent(in) :: dat !< tabulated L0L1 values
    real*8, intent(in)    :: ma     !< test particle mass [kg]
    integer*1, intent(in) :: qa     !< test particle charge number [1]
    real*8, intent(in)    :: ne     !< background electron density [1/m^3]
    real*8, intent(in)    :: the    !< normalized electron temperatures [T_b/(m_b*c^2)]
    real*8, intent(in)    :: ni(:)  !< list of background ion densities [1/m^3]
    real*8, intent(in)    :: thi(:) !< normalized ion temperatures [T_b/(m_b*c^2)]
    real*8, intent(in)    :: dt     !< time step length [s]
    real*8, intent(in)    :: rnd(3) !< array with three elements of standard normal random numbers ~ N(0,1)
    real*8, intent(in)    :: uin(3) !< normalized test particle momentum [p/mc]
    
    real*8, intent(out) :: uout(3) !< updated momentum [p/mc]

    real*8 :: clogae
    real*8, allocatable :: clogai(:)
    real*8 :: K,Dpar,Dperp,Kb,Dparb,Dperpb ! the fokker-planck coefficients
    real*8 :: dW(3)                        ! the change in the Wiener process during dt
    real*8 :: uhat(3)                      ! a unit vector parallel to pin
    real*8 :: u                            ! absolute value of particle momentum normalized to mc
    integer  :: i                          ! for iterating over plasma species

    ! Wiener process for this step
    dW = sqrt(dt)*rnd

     ! Evaluate and sum Fokker-Planck coefficients
    u = norm2(uin)
    uhat = uin / u

    allocate(clogai(size(dat%mi)))
    call ccoll_clog(ma,qa*EL_CHG,dat%mi,dat%Z0*EL_CHG,ne,the,ni,thi,u,clogae,clogai)

    ! Electron contribution
    call ccoll_coeffs(dat,ma,qa*EL_CHG,clogae,MASS_ELECTRON,-EL_CHG,ne,the,u,K=Kb,Dpar=Dparb,Dperp=Dperpb)
    K     = Kb
    Dpar  = Dparb
    Dperp = Dperpb

    ! Ion contribution
    do i = 1,size(dat%mi)
       call ccoll_coeffs(dat,ma,qa*EL_CHG,clogai(i),dat%mi(i),dat%Z0(i)*EL_CHG,ni(i),thi(i),u,K=Kb,Dpar=Dparb,Dperp=Dperpb)
       K     = K     + Kb
       Dpar  = Dpar  + Dparb
       Dperp = Dperp + Dperpb
    end do
    deallocate(clogai)

    ! Use Euler-Maruyama method to get uout
    uout = uin + K * uhat * dt + sqrt( 2.d0 * Dpar ) * dot_product( uhat, dW ) * uhat &
         + sqrt( 2.d0 * Dperp ) * ( dW - dot_product( uhat, dW ) * uhat )

  end subroutine ccoll_kinetic_relativistic_explicitpush

  
  !> Updates guiding center momentum after collisions
  !> Pushing is done by calling the explicit push function. This function is just a wrapper
  !> that additionally evaluates the plasma quantities and takes care of the coordinate transformation
  !> in momentum space thus simplifying the process of including collisions in simulations.
  subroutine ccoll_gc_relativistic_push(dat, prt, fields, mass, time, dt, jump_coll)
    implicit none
    class(ccoll_data), intent(in) :: dat !< Collision data
    class(particle_gc_relativistic), intent(inout) :: prt
    class(fields_base), intent(in) :: fields
    real*8,intent(in) :: mass, time, dt !< Mass in AMU and time in seconds
    real*8,optional,intent(in) :: jump_coll !< if provided, check if a marker perform a non-physical jump from p<p_th to p>jump_coll * p_th (used in hot_tail_tracker.f90)

    real*8 :: pnorm, E(3), B(3), psi, U, Te, Ti, the, ne, rnd(2), pin, xiin, pout, xiout, pth
    real*8, allocatable :: ni(:), thi(:)
    integer :: ierr, i_try_jump
    logical :: is_jump

    allocate(ni(size(dat%mi)), thi(size(dat%mi)))
    call fields%calc_EBpsiU(time, prt%i_elm, prt%st, prt%x(3), E, B, psi, U)                                                                                                                                                                                                          
    pnorm = sqrt(prt%p(2) * 2 * norm2(B) * mass + prt%p(1)**2)
    pin   = pnorm / ( mass * SPEED_OF_LIGHT )
    xiin  = prt%p(1) / pnorm
    call fields%calc_NjTj(time, prt%i_elm, prt%st, prt%x(3), dat%m_i_over_m_imp, ne, Te, ni, Ti)
    
    if (present(jump_coll)) then 
      ! Limit low temperature to avoid the divergence of the collision operator
      if(Te*K_BOLTZ/EL_CHG .lt. 10.0) then
         Te = 10.0*EL_CHG/K_BOLTZ
      end if
      if(Ti*K_BOLTZ/EL_CHG .lt. 10.0) then
         Ti = 10.0*EL_CHG/K_BOLTZ
      end if
      the = Te * K_BOLTZ / ( MASS_ELECTRON * SPEED_OF_LIGHT**2 )
      thi = Ti * K_BOLTZ / ( dat%mi * SPEED_OF_LIGHT**2 )

      pth = SQRT((1 + the)**2 - 1)
      is_jump = .TRUE.
      i_try_jump = 0

      ! Try the collision push and assess if a jump occured
      do while (is_jump)
         call random_number(rnd)
         rnd = floor(2.d0*rnd)
         rnd = -1.d0 + 2.d0 * rnd

         call ccoll_gc_relativistic_explicitpush(dat, mass * ATOMIC_MASS_UNIT, prt%q,  &
               ABS(ne), the, ABS(ni), thi, pin, pout, xiin, xiout, dt, rnd, 0.5*SQRT(2*the), ierr)    !<- change for different adapative reflection condition

         ! Condition of what is considered as a jump (in ||p|| or Xi)
         is_jump = ((pin .lt. pth) .AND. (pout .gt. jump_coll*pth)) .OR. (ierr .ne. 0) .OR. ((ieee_is_nan(xiout)))
         ! If not stable after 3 tries, set the out momentum to the thermal momentum isotropically
         if ((i_try_jump .gt. 3)) then
            pout = pth
            call random_number(rnd)
            xiout = 2.d0 * rnd(1) - 1.d0
            is_jump = .FALSE.
            ierr = 0
         end if 
         i_try_jump = i_try_jump + 1
      end do
    else
      the = Te * K_BOLTZ / ( MASS_ELECTRON * SPEED_OF_LIGHT**2 )
      thi = Ti * K_BOLTZ / ( dat%mi * SPEED_OF_LIGHT**2 )

      call random_number(rnd)
      rnd = floor(2.d0*rnd)
      rnd = -1.d0 + 2.d0 * rnd
      call ccoll_gc_relativistic_explicitpush(dat, mass * ATOMIC_MASS_UNIT, prt%q,  &
            ne, the, ni, thi, pin, pout, xiin, xiout, dt, rnd, 1.0e-4, ierr)
    end if

    if (ierr .ne. 0) then
      prt%i_elm = 0
      write(*,*) "WARNING: |xiout| > 2 after particle collision, removing particle. Try reducing particle time step"
    end if

    pnorm = pout * ( mass * SPEED_OF_LIGHT )
    prt%p(1) = pnorm * xiout
    prt%p(2) = ( pnorm**2 - prt%p(1)**2 ) / ( 2 * norm2(B) * mass )
    deallocate(ni, thi)

  end subroutine ccoll_gc_relativistic_push

  ! Apply Coulomb collisions for guiding center
  subroutine ccoll_gc_relativistic_explicitpush(dat,ma,qa,ne,the,ni,thi,uin,uout,xiin,xiout,dt,rnd,cutoff,ierr)
    
    class(ccoll_data), intent(in) :: dat !< tabulated L0L1 values
    real*8, intent(in)    :: ma     !< test particle mass [kg]
    integer*1, intent(in) :: qa     !< test particle charge number [1]
    real*8, intent(in)    :: ne     !< background electron density [1/m^3]
    real*8, intent(in)    :: the    !< normalized electron temperature [T_b/(m_b*c^2)]
    real*8, intent(in)    :: ni(:)  !< list of background ion densities [1/m^3]
    real*8, intent(in)    :: thi(:) !< normalized ion temperatures [T_b/(m_b*c^2)]
    real*8, intent(in)    :: dt     !< time step length [s]
    real*8, intent(in)    :: uin    !< test particle momentum  [p/mc]
    real*8, intent(in)    :: xiin   !< test particle pitch [ppar/p]
    real*8, intent(in)    :: cutoff !< minimum normalized momentum, energies below this are reflected
    real*8, intent(in)    :: rnd(2) !< normally ditributed random numbes

    real*8,   intent(out) :: uout  !< updated momentum
    real*8,   intent(out) :: xiout !< updated pitch
    integer,  intent(out) :: ierr

    real*8 :: clogae
    real*8, allocatable :: clogai(:)
    real*8 :: kappa, Dpar, dDpar, Dperp, nu ! Collision coefficients
    real*8 :: kappab, Dparb, dDparb, Dperpb ! Coll. coefficients species-wise
    integer :: i

    ierr = 0

    allocate(clogai(size(dat%mi)))
    call ccoll_clog(ma,qa*EL_CHG,dat%mi,dat%Z0*EL_CHG,ne,the,ni,thi,uin,clogae,clogai)

    ! Electron contribution
    call ccoll_coeffs(dat, ma, qa*EL_CHG, clogae, MASS_ELECTRON, -EL_CHG, ne, the, uin, &
            kappa=kappab, Dpar=Dparb, Dperp=Dperpb, dDpar=dDparb)
    kappa = kappab
    Dpar  = Dparb
    dDpar = dDparb
    Dperp = Dperpb

    ! Ion contribution
    do i = 1,size(dat%mi)
       call ccoll_coeffs(dat, ma, qa*EL_CHG, clogai(i), dat%mi(i), dat%Z0(i)*EL_CHG, ni(i), thi(i), uin, &
            kappa=kappab, Dpar=Dparb, Dperp=Dperpb, dDpar=dDparb)

       kappa = kappa + kappab
       Dpar  = Dpar  + Dparb
       dDpar = dDpar + dDparb
       Dperp = Dperp + Dperpb
    end do

    nu = 2.d0 * Dperp / uin**2
    deallocate(clogai)

    uout  = uin + ( kappa + dDpar + 2.d0 * Dpar / uin ) * dt &
                + sqrt( 2.d0 * Dpar * dt ) * rnd(1)
    xiout = xiin - xiin * nu * dt + sqrt( ( 1.0 - xiin**2 ) * nu * dt ) * rnd(2)

    ! Reflect uout if uout is below the cutoff value
    if(uout .lt. cutoff) then
       uout = 2.d0 * cutoff-uout
    end if

    ! Reflect pitch if xiout is outside the interval [-1, 1]
    if(abs(xiout) .gt. 1.D0) then
       ! First make sure xiout is between the interval [-2, 2]. Physics-wise what we do here is not justified,
       ! but neither is having |xiout| > 2 (one should decrease time step if this happens).
       if (abs(xiout) .gt. 2.D0) ierr = 1

       xiout = mod( xiout, 2.0 )

       ! Reflect (this part is ok physics-wise)
       if(abs(xiout) .gt. 1.D0) then
          xiout = sign(2.D0-abs(xiout), xiout)
       end if
    end if

  end subroutine ccoll_gc_relativistic_explicitpush

  
  !> Evaluate coefficients for the partial screening operator
  !> Note that the partial screening operator assumes test particle is electron.
  subroutine ccoll_coeffs_partialscreening(dat, ne, the, ni, u, nuee_S, nuee_par, nuee_D, nuei_D)
    implicit none

    type(ccoll_data), intent(in) :: dat
    real*8, intent(in)  :: ne
    real*8, intent(in)  :: the
    real*8, intent(in)  :: ni(:)
    real*8, intent(in)  :: u
    real*8, intent(out) :: nuee_S
    real*8, intent(out) :: nuee_par
    real*8, intent(out) :: nuee_D
    real*8, intent(out) :: nuei_D

    real*8 :: clog0, clogee, clogei
    real*8 :: gamma, nu_c, Zeff
    real*8 :: mu0, mu1, mu2
    real*8 :: Nebnd, g, h, tmp
    integer :: i

    gamma = sqrt(1.d0 + u**2)
    call ccoll_mufuncs(dat,u,the,mu0,mu1,mu2)

    ! Coulomb logarithms specifically for energetic electrons
    clog0 = 14.9d0 - 0.5d0 * log( ne / 1.d20 ) + log( the * MASS_ELECTRON * SPEED_OF_LIGHT**2 / ( 1.d3 * EL_CHG ) ) 
    clogee = clog0 + (1.d0/5.d0) * log( 1.d0 + ( 2.d0 * ( gamma - 1.d0 ) / sqrt( 2.d0 * the )  )**(5.d0/2.d0) )
    clogei = clog0 + (1.d0/5.d0) * log( 1.d0 + ( 2.d0 * u / sqrt( 2.d0 * the ) )**5 )

    Zeff = 0.d0
    g    = 0.d0
    h    = 0.d0
    do i=1,size(ni)
       Nebnd = dat%Zi(i) - dat%Z0(i)
       Zeff  = Zeff + ni(i) * dat%Z0(i)**2 / ne

       tmp = ( u * dat%ai(i) )**( 3.d0 / 2.d0 )
       g = g + ( 2.d0 / 3.d0 ) * ( ni(i) / ne ) * ( ( dat%Zi(i)**2 - dat%Z0(i)**2 ) * log( tmp + 1.d0 ) &
             - Nebnd**2 * ( tmp / ( tmp + 1.d0 ) ) )

       h = h + ( ni(i) / ne ) * Nebnd * ( log( 1.d0 + ( u * sqrt( gamma - 1.d0 ) / ( dat%Ii(i) / ( MASS_ELECTRON * SPEED_OF_LIGHT**2 ) ) )**5 ) / 5.d0 &
             - 1.d0 + 1.d0 / gamma**2 )
    end do
    
    nu_c     = ( ne * EL_CHG**4 * clog0 ) / ( 4 * PI * EPS_ZERO**2 * MASS_ELECTRON**2 * SPEED_OF_LIGHT**3 )
    nuee_S   = ( nu_c * gamma**2 / ( clog0 * u**2 ) ) * ( clogee * mu1 / gamma**2 + h )
    nuee_par = ( 2 * nu_c * gamma * the / u**3 ) * mu1
    nuee_D   = ( 2 * nu_c * clogee / ( u**2 * clog0 ) ) * ( u**2 * mu0 + u**2 * gamma * the * mu2 - the * mu1 ) / ( 2 * gamma * u**3 )
    nuei_D   = ( nu_c * gamma / ( u**3 * clog0 ) ) * ( clogei * Zeff + g )

  end subroutine ccoll_coeffs_partialscreening

  !> Push electron taking partial screening effect into account
  subroutine ccoll_explicitpush_partialscreening(dat, ne, the, ni, uin, uout, xiin, xiout, dt, rnd, cutoff, ierr)
    implicit none
    type(ccoll_data), intent(in) :: dat
    
    real*8, intent(in) :: ne, ni(:),the, uin, xiin
    real*8, intent(in) :: cutoff, dt, rnd(2)
    real*8, intent(out) :: uout, xiout
    integer, intent(out) :: ierr

    real*8 :: nuee_S, nuee_par, nuee_D, nuei_D

    ierr = 0

    call ccoll_coeffs_partialscreening(dat, ne, the, ni, uin, nuee_S, nuee_par, nuee_D, nuei_D)

    
    uout  = uin - nuee_S * dt + sqrt( nuee_par * dt ) * rnd(1)
    if(uout .lt. cutoff) then
       uout = 2*cutoff-uout
    end if
    
    xiout = xiin - ( nuee_D + nuei_D ) * xiin * dt + sqrt( ( 1.d0 - xiin**2 ) * ( nuee_D + nuei_D ) * dt  ) * rnd(2)
    if(abs(xiout) .gt. 1.d0) then
       if (abs(xiout) .gt. 2.D0) ierr = 1
       xiout = mod( xiout, 2.d0 )
       if(abs(xiout) .gt. 1.d0) then
          xiout = sign(2.d0-abs(xiout), xiout)
       end if
    end if

  end subroutine ccoll_explicitpush_partialscreening

  !> Push guiding center electron taking partial screening into account
  !> Evaluates the field and takes care of the coordinate transformation before calling the
  !> explicit pusher.
  subroutine ccoll_gc_relativistic_push_partialscreening(dat, prt, fields, mass, time, dt)
    implicit none
    class(ccoll_data), intent(in) :: dat !< Collision data
    class(particle_gc_relativistic), intent(inout) :: prt
    class(fields_base), intent(in) :: fields
    real*8,intent(in) :: mass, time, dt !< Mass in AMU and time in seconds
    real*8 :: pnorm, E(3), B(3), psi, U, Te, Ti, the, ne, rnd(2), pin, xiin, pout, xiout
    real*8, allocatable :: ni(:), thi(:)
    integer :: ierr

    allocate(ni(size(dat%mi)), thi(size(dat%mi)))
    call fields%calc_EBpsiU(time, prt%i_elm, prt%st, prt%x(3), E, B, psi, U)
    pnorm = sqrt(prt%p(2) * 2 * norm2(B) * mass + prt%p(1)**2)
    pin   = pnorm / ( mass * SPEED_OF_LIGHT )
    xiin  = prt%p(1) / pnorm
    call fields%calc_NjTj(time, prt%i_elm, prt%st, prt%x(3), dat%m_i_over_m_imp, ne, Te, ni, Ti)
    the = Te * K_BOLTZ / ( MASS_ELECTRON * SPEED_OF_LIGHT**2 )
    thi = Ti * K_BOLTZ / ( dat%mi * SPEED_OF_LIGHT**2 )

    call random_number(rnd)
    rnd = floor(2.d0*rnd)
    rnd = -1.d0 + 2.d0 * rnd

    call ccoll_explicitpush_partialscreening(dat, ne, the, ni, pin, pout, xiin, xiout, dt, rnd, DEFAULT_UCUTOFF, ierr)
     
    if (ierr .ne. 0) then
      prt%i_elm = 0
      write(*,*) "WARNING: |xiout| > 2 after particle collision, removing particle. Try reducing particle time step"
    end if

    pnorm = pout * ( mass * SPEED_OF_LIGHT )
    prt%p(1) = pnorm * xiout
    prt%p(2) = ( pnorm**2 - prt%p(1)**2 ) / ( 2 * norm2(B) * mass )
    deallocate(ni, thi)

  end subroutine ccoll_gc_relativistic_push_partialscreening

  !> Push gyro orbiting electron taking partial screening into account
  !> Evaluates the field and takes care of the coordinate transformation before calling the
  !> explicit pusher.
  subroutine ccoll_kinetic_relativistic_push_partialscreening(dat, prt, fields, mass, time, dt)
    implicit none
    class(ccoll_data), intent(in) :: dat !< Collision data
    class(particle_kinetic_relativistic), intent(inout) :: prt
    class(fields_base), intent(in) :: fields
    real*8,intent(in) :: mass, time, dt !< Mass in AMU and time in seconds
    real*8 :: E(3), B(3), psi, U, ne, rnd(2), pin, pout, xiin, xiout, Te, Ti, the, bperp(3), bhat(3)
    real*8, allocatable :: ni(:), thi(:)
    integer :: ierr

    call fields%calc_EBpsiU(time, prt%i_elm, prt%st, prt%x(3), E, B, psi, U)
    bhat = vector_cylindrical_to_cartesian(prt%x(3), B) / norm2(B)
    pin  = norm2(prt%p) / (mass * SPEED_OF_LIGHT)
    xiin = dot_product(prt%p, bhat) / ( norm2(prt%p) )

    allocate(ni(size(dat%mi)), thi(size(dat%mi)))
    call fields%calc_NjTj(time, prt%i_elm, prt%st, prt%x(3), dat%m_i_over_m_imp, ne, Te, ni, Ti)
    the = Te * K_BOLTZ / ( MASS_ELECTRON * SPEED_OF_LIGHT**2 )
    thi = Ti * K_BOLTZ / ( dat%mi * SPEED_OF_LIGHT**2 )

    call random_number(rnd)
    rnd = floor(2.d0*rnd)
    rnd = -1.d0 + 2.d0 * rnd

    call ccoll_explicitpush_partialscreening(dat, ne, the, ni, pin, pout, xiin, xiout, dt, rnd, DEFAULT_UCUTOFF, ierr)

    if (ierr .ne. 0) then
      prt%i_elm = 0
      write(*,*) "WARNING: |xiout| > 2 after particle collision, removing particle. Try reducing particle time step"
    end if
    
    ! Back to particle coordinates (gyroangle is left invariant)
    bperp = prt%p - dot_product(prt%p, bhat) * bhat
    bperp = bperp / norm2(bperp)

    prt%p = (xiout * bhat + sqrt( 1.d0 - xiout**2 ) * bperp ) * pout * (mass * SPEED_OF_LIGHT)
    deallocate(ni, thi)

  end subroutine ccoll_kinetic_relativistic_push_partialscreening

  !> Evaluates the mu functions (and their derivatives if needed)
  subroutine ccoll_mufuncs(data,u,th,mu0,mu1,mu2,dmu0,dmu1,dmu2)
    implicit none
    class(ccoll_data), intent(in) :: data !< initialized L0L1 tables
    real*8, intent(in)  :: u   !< p/mc value
    real*8, intent(in)  :: th  !< T/mc^2
    real*8, intent(out) :: mu0 !< Eq. 14 in the referece paper
    real*8, intent(out) :: mu1 !< Eq. 15 in the referece paper
    real*8, intent(out) :: mu2 !< Eq. 16 in the referece paper
    real*8, intent(out), optional :: dmu0 !< d mu0 / d u
    real*8, intent(out), optional :: dmu1 !< d mu1 / d u
    real*8, intent(out), optional :: dmu2 !< d mu2 / d u

    real*8 :: gamma,gammasq,expBessel2,L0,L1,expgammatheta,tg,th2,u2,tgK

    u2            = u**2
    gammasq       = 1.D0 + u2
    gamma         = sqrt(gammasq)
    expBessel2    = bessel_k2exp(1.D0 / th)
    expgammatheta = exp( ( 1.D0 - gamma ) / th )
    tg  = th * gamma
    th2 = th**2

    call interp_L0L1(data,u,th,L0,L1)

    mu0 = ( gammasq * L0 - th * L1 + ( th - gamma ) * u * expgammatheta ) / expBessel2
    mu1 = ( gammasq * L1 - th * L0 + ( th * gamma - 1.d0 ) * u * expgammatheta ) / expBessel2
    mu2 = ( 2.d0 * tg * L1 + ( 1 + 2.d0 * th2 ) * u * expgammatheta ) / ( th * expBessel2 )

    if(present(dmu0) .or. present(dmu1) .or. present(dmu2)) then
       tgK = tg * expBessel2
       if(present(dmu0)) then
          dmu0 = ( 2.d0 * tg * u * L0 + ( gamma - 2.d0 * th ) * u2 * expgammatheta ) / tgK
       end if
       if(present(dmu1)) then
          !dmu1=(2.D0*tg*u*L1+(2.D0*th2+1.D0)*u2*expgammatheta)/tgK ! This is the explicit form
          dmu1 = mu2 * u / gamma
       end if
       if(present(dmu2)) then
          dmu2 = ( 2.d0 * th2 * u * L1 + ( 2.d0 * th2 * tg + 2.d0 * th2 + tg - u2 ) * expgammatheta ) / ( tgK * th )
       end if
    end if
    
  end subroutine ccoll_mufuncs


  !> Interpolates (bilinear) the special functions L0 and L1
  !> Approximations are used outside the tabulated domain.
  subroutine interp_L0L1(data,u,theta,L0,L1)
    class(ccoll_data), intent(in) :: data !< tabulated L0L1 values
    real*8, intent(in)  :: u     !< queried p/mc value
    real*8, intent(in)  :: theta !< queried T/mc^2 value
    real*8, intent(out) :: L0    !< interpolated L0
    real*8, intent(out) :: L1    !< interpolated L1

    real*8  :: th
    integer :: nu, nth

    nu = size(data%u)
    nth = size(data%theta)

    if ( theta .gt. data%theta(nth) ) then
       th=data%theta(nth) 
       write(*,*) 'Warning: temperature exceeds tabulated', theta, th
    else
       th=theta
    end if
    
    if(u.gt.data%u(nu)) then
       L0 = bessel_k0exp(1.D0/th)
       L1 = bessel_k1exp(1.D0/th)
    else
       if( (u.gt.data%u(1)) .and. (th.gt.data%theta(1)) ) then
          
          L0 = interp_bilinear(log10(data%u),log10(data%theta),data%L0,log10(u),log10(th)) 
          L1 = interp_bilinear(log10(data%u),log10(data%theta),data%L1,log10(u),log10(th))
       else
          L0 = sqrt(pi*theta/2)*erf(u/sqrt(2*theta))
          L1 = L0
       end if
    end if
    
  end subroutine interp_L0L1


end module mod_ccoll_relativistic


