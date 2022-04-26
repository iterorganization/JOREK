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

xout(1:n)  =  (1.d0 - s)**2  * (      (1.d0 - s) * xx(1,1:n)  + 3.d0 * s  * xx(2,1:n) ) &
           +   s**2          * ( 3.d0*(1.d0 - s) * xx(3,1:n)  +        s  * xx(4,1:n) )

return
end
