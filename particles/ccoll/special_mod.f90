!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> @brief
!> A module for evaluating special functions needed in 
!> computing the relativistic collision coefficients.
!>
!> contains:
!> \f$ L_0(x,\Theta)=\int_0^{x}ds\frac{e^{(1-\sqrt{1+s^2})/\Theta}}{\sqrt{1+s^2}} \f$ \\
!> \f$ L_1(x,\Theta)=\int_0^{x}ds\; e^{(1-\sqrt{1+s^2})/\Theta} \f$ \\
!> \f$ exp(x)K_2(x) \f$
!<
module special_mod 
  use func_mod, only : func_real8_1d
  use simpson_mod, only : simpson_adaptive
  implicit none
  
  private

  real*8, parameter :: DEFAULT_L0L1_eps    = 1.D-8
  real*8, parameter :: DEFAULT_L0L1_cutoff = 1.D-7

  public :: special_L0L1, &
           &special_besk2exp, &
           &special_besk1exp, &
           &special_besk0exp, &
           &special_test

contains

  !> Evaluates the functions \f$ L_0(x,\Theta), L_1(x,\Theta) \f$0
  !> The cutoff parameter (< 1) divides the integral in two
  !> parts at x_cutoff = sqrt((1-theta * ln(cutoff))**2 -1).
  !> This helps to ensure that adaptive integration does not
  !> fail when theta is very small. 
  subroutine special_L0L1(x,theta,L0,L1,L0L1_eps,L0L1_cutoff)
    implicit none
    real*8, intent(in) :: x,theta
    real*8, intent(in), optional :: L0L1_eps,L0L1_cutoff

    real*8, intent(out) :: L0,L1

    real*8 :: def_cutoff, x_cutoff
    real*8 :: def_tol, tol
    procedure(func_real8_1d), pointer :: L0_ptr => null(), L1_ptr => null()

    def_tol = DEFAULT_L0L1_eps
    def_cutoff = DEFAULT_L0L1_cutoff
    if(present(L0L1_eps)) def_tol = L0L1_eps
    if(present(L0L1_cutoff)) def_cutoff = L0L1_cutoff

    x_cutoff=sqrt((1.D0 - theta * log(def_cutoff))**2 - 1.D0)

    L0_ptr => L0_integrand
    L1_ptr => L1_integrand
    
    if(x_cutoff > x) then
       tol = 0.5D0 * def_tol * x
       L0 = simpson_adaptive(L0_ptr,0.D0,x,tol,10)
       L1 = simpson_adaptive(L1_ptr,0.D0,x,tol,10)
    else
       tol = 0.5D0 * def_tol * x_cutoff
       L0 = simpson_adaptive(L0_ptr,0.D0,x_cutoff,tol,10)+simpson_adaptive(L0_ptr,x_cutoff,x,tol,20)
       L1 = simpson_adaptive(L1_ptr,0.D0,x_cutoff,tol,10)+simpson_adaptive(L1_ptr,x_cutoff,x,tol,20)
    end if

  contains

    function L0_integrand(x) result (val)
      real*8, intent(in) :: x
      real*8 :: val,gamma
      gamma = sqrt(1.D0+x**2)
      val = exp((1.D0-gamma)/theta)/gamma
    end function L0_integrand
    
    function L1_integrand(x) result (val)
      real*8, intent(in) :: x
      real*8 :: val,gamma
      gamma = sqrt(1.D0+x**2)
      val = exp((1.D0-gamma)/theta)
    end function L1_integrand

  end subroutine special_L0L1

  !> Few functions for evaluating modified Bessel function of second type
  !> Source: http://people.sc.fsu.edu/~jburkardt/f_src/specfun/specfun.f90
  !> besk2 evaluates function K2. Higher order functions could be evaluated
  !> easily by implementing this function recursively.
  !> 
  !> NOTE: output of functions special_besk0, special_besk1,
  !> and special_besk2 is scaled
  !> by scalar e^(x), i.e., f_scaled = e^(x)*f_exact 
  !> to avoid under & overflow
  !> issues.
  !<
  function special_besk2exp(x) result(out)

    implicit none

    real (kind = 8) x,out

    if(x .eq. 0.0) then
       out = special_besk0exp(x)
    else
       out = special_besk0exp(x)+2*special_besk1exp(x)/x 
    end if

  end function special_besk2exp

  !*****************************************************************************80
  !
  !! BESK0 evaluates the Bessel K0(X) function.
  !
  !  Discussion:
  !
  !    This routine computes approximate values for the
  !    modified Bessel function of the second kind of order zero
  !    for arguments 0.0 < ARG <= XMAX.
  !
  !    See comments heading CALCK0.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 April 2007
  !
  !  Author:
  !
  !    Original FORTRAN77 version by William Cody, Laura Stoltz.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) X, the argument of the function.
  !
  !    Output, real ( kind = 8 ) BESK0, the value of the function.
  !

  function special_besk0exp ( x )

    implicit none

    real ( kind = 8 ) special_besk0exp
    integer ( kind = 4 ) jint
    real ( kind = 8 ) result
    real ( kind = 8 ) x

    jint = 2 ! Scaled = true
    call special_calck0 ( x, result, jint )
    special_besk0exp = result

    return
  end function special_besk0exp
  

  !*****************************************************************************80
  !
  !! BESK1 evaluates the Bessel K1(X) function.
  !
  !  Discussion:
  !
  !    This routine computes approximate values for the
  !    modified Bessel function of the second kind of order one
  !    for arguments XLEAST <= ARG <= XMAX.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    03 April 2007
  !
  !  Author:
  !
  !    Original FORTRAN77 version by William Cody.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) X, the argument of the function.
  !
  !    Output, real ( kind = 8 ) BESK1, the value of the function.
  !

  function special_besk1exp ( x )

    implicit none

    real ( kind = 8 ) special_besk1exp
    integer ( kind = 4 ) jint
    real ( kind = 8 ) result
    real ( kind = 8 ) x

    jint = 2 ! Scaled = true
    call special_calck1 ( x, result, jint )
    special_besk1exp = result

    return
  end function special_besk1exp


  subroutine special_calck0 ( arg, result, jint )

    !*****************************************************************************80
    !
    !! CALCK0 computes various K0 Bessel functions.
    !
    !  Discussion:
    !
    !    This routine computes modified Bessel functions of the second kind
    !    and order zero, K0(X) and EXP(X)*K0(X), for real
    !    arguments X.
    !
    !    The main computation evaluates slightly modified forms of near
    !    minimax rational approximations generated by Russon and Blair,
    !    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
    !    1969.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) ARG, the argument.  0 < ARG is
    !    always required.  If JINT = 1, then the argument must also be
    !    less than XMAX.
    !
    !    Output, real ( kind = 8 ) RESULT, the value of the function,
    !    which depends on the input value of JINT:
    !    1, RESULT = K0(x);
    !    2, RESULT = exp(x) * K0(x);
    !
    !    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
    !    1, K0(x);
    !    2, exp(x) * K0(x);
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) jint
    real ( kind = 8 ) arg
    real ( kind = 8 ) f(4)
    real ( kind = 8 ) g(3)
    real ( kind = 8 ) p(6)
    real ( kind = 8 ) pp(10)
    real ( kind = 8 ) q(2)
    real ( kind = 8 ) qq(10)
    real ( kind = 8 ) result
    real ( kind = 8 ) sumf
    real ( kind = 8 ) sumg
    real ( kind = 8 ) sump
    real ( kind = 8 ) sumq
    real ( kind = 8 ) temp
    real ( kind = 8 ) x
    real ( kind = 8 ) xinf
    real ( kind = 8 ) xmax
    real ( kind = 8 ) xsmall
    real ( kind = 8 ) xx
    !
    !  Machine-dependent constants
    !
    data xsmall /1.11d-16/
    data xinf /1.79d+308/
    data xmax /705.342d0/
    !
    !  Coefficients for XSMALL <= ARG <= 1.0
    !
    data   p/ 5.8599221412826100000d-04, 1.3166052564989571850d-01, &
         1.1999463724910714109d+01, 4.6850901201934832188d+02, &
         5.9169059852270512312d+03, 2.4708152720399552679d+03/
    data   q/-2.4994418972832303646d+02, 2.1312714303849120380d+04/
    data   f/-1.6414452837299064100d+00,-2.9601657892958843866d+02, &
         -1.7733784684952985886d+04,-4.0320340761145482298d+05/
    data   g/-2.5064972445877992730d+02, 2.9865713163054025489d+04, &
         -1.6128136304458193998d+06/
    !
    !  Coefficients for  1.0 < ARG
    !
    data  pp/ 1.1394980557384778174d+02, 3.6832589957340267940d+03, &
         3.1075408980684392399d+04, 1.0577068948034021957d+05, &
         1.7398867902565686251d+05, 1.5097646353289914539d+05, &
         7.1557062783764037541d+04, 1.8321525870183537725d+04, &
         2.3444738764199315021d+03, 1.1600249425076035558d+02/
    data  qq/ 2.0013443064949242491d+02, 4.4329628889746408858d+03, &
         3.1474655750295278825d+04, 9.7418829762268075784d+04, &
         1.5144644673520157801d+05, 1.2689839587977598727d+05, &
         5.8824616785857027752d+04, 1.4847228371802360957d+04, &
         1.8821890840982713696d+03, 9.2556599177304839811d+01/

    x = arg
    !
    !  0.0 < ARG <= 1.0.
    !
    if ( 0.0D+00 < x ) then

       if ( x <= 1.0D+00 ) then

          temp = log ( x )

          if ( x < xsmall ) then
             !
             !  Return for small ARG.
             !
             result = p(6) / q(2) - temp

          else

             xx = x * x

             sump = (((( &
                  p(1) &
                  * xx + p(2) ) &
                  * xx + p(3) ) &
                  * xx + p(4) ) &
                  * xx + p(5) ) &
                  * xx + p(6)

             sumq = ( xx + q(1) ) * xx + q(2)
             sumf = ( ( &
                  f(1) &
                  * xx + f(2) ) &
                  * xx + f(3) ) &
                  * xx + f(4)

             sumg = ( ( xx + g(1) ) * xx + g(2) ) * xx + g(3)

             result = sump / sumq - xx * sumf * temp / sumg - temp

             if ( jint == 2 ) then
                result = result * exp ( x )
             end if

          end if

       else if ( jint == 1 .and. xmax < x ) then
          !
          !  Error return for XMAX < ARG.
          !
          result = 0.0D+00

       else
          !
          !  1.0 < ARG.
          !
          xx = 1.0D+00 / x
          sump = pp(1)
          do i = 2, 10
             sump = sump * xx + pp(i)
          end do

          sumq = xx
          do i = 1, 9
             sumq = ( sumq + qq(i) ) * xx
          end do
          sumq = sumq + qq(10)
          result = sump / sumq / sqrt ( x )

          if ( jint == 1 ) then
             result = result * exp ( -x )
          end if

       end if

    else
       !
       !  Error return for ARG <= 0.0.
       !
       result = xinf

    end if

    return
  end subroutine special_calck0

  subroutine special_calck1 ( arg, result, jint )

    !*****************************************************************************80
    !
    !! CALCK1 computes various K1 Bessel functions.
    !
    !  Discussion:
    !
    !    This routine computes modified Bessel functions of the second kind
    !    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
    !
    !    The main computation evaluates slightly modified forms of near
    !    minimax rational approximations generated by Russon and Blair,
    !    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
    !    1969.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) ARG, the argument.  XLEAST < ARG is
    !    always required.  If JINT = 1, then the argument must also be
    !    less than XMAX.
    !
    !    Output, real ( kind = 8 ) RESULT, the value of the function,
    !    which depends on the input value of JINT:
    !    1, RESULT = K1(x);
    !    2, RESULT = exp(x) * K1(x);
    !
    !    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
    !    1, K1(x);
    !    2, exp(x) * K1(x);
    !
    implicit none

    real ( kind = 8 ) arg
    real ( kind = 8 ) f(5)
    real ( kind = 8 ) g(3)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) jint
    real ( kind = 8 ) p(5)
    real ( kind = 8 ) pp(11)
    real ( kind = 8 ) q(3)
    real ( kind = 8 ) qq(9)
    real ( kind = 8 ) result
    real ( kind = 8 ) sumf
    real ( kind = 8 ) sumg
    real ( kind = 8 ) sump
    real ( kind = 8 ) sumq
    real ( kind = 8 ) x
    real ( kind = 8 ) xinf
    real ( kind = 8 ) xmax
    real ( kind = 8 ) xleast
    real ( kind = 8 ) xsmall
    real ( kind = 8 ) xx
    !
    !  Machine-dependent constants
    !
    data xleast /2.23d-308/
    data xsmall /1.11d-16/
    data xinf /1.79d+308/
    data xmax /705.343d+0/
    !
    !  Coefficients for  XLEAST <=  ARG  <= 1.0
    !
    data   p/ 4.8127070456878442310d-1, 9.9991373567429309922d+1, &
         7.1885382604084798576d+3, 1.7733324035147015630d+5, &
         7.1938920065420586101d+5/
    data   q/-2.8143915754538725829d+2, 3.7264298672067697862d+4, &
         -2.2149374878243304548d+6/
    data   f/-2.2795590826955002390d-1,-5.3103913335180275253d+1, &
         -4.5051623763436087023d+3,-1.4758069205414222471d+5, &
         -1.3531161492785421328d+6/
    data   g/-3.0507151578787595807d+2, 4.3117653211351080007d+4, &
         -2.7062322985570842656d+6/
    !
    !  Coefficients for  1.0 < ARG
    !
    data  pp/ 6.4257745859173138767d-2, 7.5584584631176030810d+0, &
         1.3182609918569941308d+2, 8.1094256146537402173d+2, &
         2.3123742209168871550d+3, 3.4540675585544584407d+3, &
         2.8590657697910288226d+3, 1.3319486433183221990d+3, &
         3.4122953486801312910d+2, 4.4137176114230414036d+1, &
         2.2196792496874548962d+0/
    data  qq/ 3.6001069306861518855d+1, 3.3031020088765390854d+2, &
         1.2082692316002348638d+3, 2.1181000487171943810d+3, &
         1.9448440788918006154d+3, 9.6929165726802648634d+2, &
         2.5951223655579051357d+2, 3.4552228452758912848d+1, &
         1.7710478032601086579d+0/

    x = arg
    !
    !  Error return for ARG < XLEAST.
    !
    if ( x < xleast ) then

       result = xinf
       !
       !  XLEAST <= ARG <= 1.0.
       !
    else if ( x <= 1.0D+00 ) then

       if ( x < xsmall ) then
          !
          !  Return for small ARG.
          !
          result = 1.0D+00 / x

       else

          xx = x * x

          sump = (((( &
               p(1) &
               * xx + p(2) ) &
               * xx + p(3) ) &
               * xx + p(4) ) &
               * xx + p(5) ) &
               * xx + q(3)

          sumq = (( &
               xx + q(1) ) &
               * xx + q(2) ) &
               * xx + q(3)

          sumf = ((( &
               f(1) &
               * xx + f(2) ) &
               * xx + f(3) ) &
               * xx + f(4) ) &
               * xx + f(5)

          sumg = (( &
               xx + g(1) ) &
               * xx + g(2) ) &
               * xx + g(3)

          result = ( xx * log ( x ) * sumf / sumg + sump / sumq ) / x

          if ( jint == 2 ) then
             result = result * exp ( x )
          end if

       end if

    else if ( jint == 1 .and. xmax < x ) then
       !
       !  Error return for XMAX < ARG.
       !
       result = 0.0D+00

    else
       !
       !  1.0 < ARG.
       !
       xx = 1.0D+00 / x

       sump = pp(1)
       do i = 2, 11
          sump = sump * xx + pp(i)
       end do

       sumq = xx
       do i = 1, 8
          sumq = ( sumq + qq(i) ) * xx
       end do
       sumq = sumq + qq(9)

       result = sump / sumq / sqrt ( x )

       if ( jint == 1 ) then
          result = result * exp ( -x )
       end if

    end if

    return
  end subroutine special_calck1



  !> a test routine that check the evaluation of the special functions
  !< 
  subroutine special_test()
    implicit none

    ! L0 and L1 functions
    real*8 :: x, theta, L0,L1,expK2
    integer:: i
    
    x=1.D1
    theta = 1.D-5
    call special_L0L1(x,theta,L0,L1)
    write(*,*) "x = ",x," theta = ",theta
    write(*,*) "(L1-L0)/L1 = ", (L1-L0)/L1
    write(*,*) "(L1-app1)/L1 = ",(L1-sqrt(2.D0*atan(1.D0)*theta)*erf(x/sqrt(2*theta)))/L1
    write(*,*) "(L1-app2)/L1 = ",(L1-special_besk1exp(1.D0/theta))/L1
    
    x=1.D-3
    theta = 1.D-5
    call special_L0L1(x,theta,L0,L1)
    write(*,*) "x = ",x," theta = ",theta
    write(*,*) "(L1-L0)/L1 = ", (L1-L0)/L1
    write(*,*) "(L1-app1)/L1 = ",(L1-sqrt(2.D0*atan(1.D0)*theta)*erf(x/sqrt(2*theta)))/L1
    write(*,*) "(L1-app2)/L1 = ",(L1-special_besk1exp(1.D0/theta))/L1

    x=1.D1
    theta = 1.D-1
    call special_L0L1(x,theta,L0,L1)
    write(*,*) "x = ",x," theta = ",theta
    write(*,*) "(L1-L0)/L1 = ", (L1-L0)/L1
    write(*,*) "(L1-app1)/L1 = ",(L1-sqrt(2.D0*atan(1.D0)*theta)*erf(x/sqrt(2*theta)))/L1
    write(*,*) "(L1-app2)/L1 = ",(L1-special_besk1exp(1.D0/theta))/L1
    write(*,*) "(L0-app3)/L0 = ",(L0-special_besk0exp(1.D0/theta))/L0
    write(*,*) "L0, exp(1/theta)BesselK(0,1/theta) = ", L0,special_besk0exp(1.D0/theta)

    x=1.D-3
    theta = 1.D-1
    call special_L0L1(x,theta,L0,L1)
    write(*,*) "x = ",x," theta = ",theta
    write(*,*) "(L1-L0)/L1 = ", (L1-L0)/L1
    write(*,*) "(L1-app1)/L1 = ",(L1-sqrt(2.D0*atan(1.D0)*theta)*erf(x/sqrt(2*theta)))/L1
    write(*,*) "(L1-app2)/L1 = ",(L1-special_besk1exp(1.D0/theta))/L1

  end subroutine special_test

end module special_mod
