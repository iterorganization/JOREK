subroutine Bezier_1d(n,s,xx,xout)
!--------------------------------------------------------------
! defined on the interval ( 0 < s < 1 )
!--------------------------------------------------------------
use mod_parameters, only: n_order
implicit none
integer :: n
real*8  :: xx(n_order+1,n)
real*8  :: xout(n),s

!xout       =             (1. - s)**3   * xx(1,:)  &
!           + 3. * s    * (1. - s)**2   * xx(2,:)  &
!           + 3. * s**2 * (1. - s)      * xx(3,:)  &
!           +      s**3                 * xx(4,:)

if (n_order .eq. 3) then
  xout(1:n)  =  (1.d0 - s)**2  * (      (1.d0 - s) * xx(1,1:n)  + 3.d0 * s  * xx(2,1:n) ) &
             +   s**2          * ( 3.d0*(1.d0 - s) * xx(3,1:n)  +        s  * xx(4,1:n) )
else
  xout(1:n)  =              (1.-s)**5 * xx(1,1:n) &
             + 5.  * s**1 * (1.-s)**4 * xx(2,1:n) &
             + 10. * s**2 * (1.-s)**3 * xx(3,1:n) &
             + 10. * s**3 * (1.-s)**2 * xx(4,1:n) &
             + 5.  * s**4 * (1.-s)**1 * xx(5,1:n) &
             +       s**5             * xx(6,1:n)
endif

return
end
