!> Contains the Bezier and Fourier basis functions on Gaussian points.
!!
!! - Array dimensions of Bezier basis functions and derivatives
!!   (H, H_s, H_t, H_st, H_ss, H_tt):
!!   - Vertex
!!   - Basis function (i.e. p0,u,v,w)
!!   - Gaussian point in s direction
!!   - Gaussian point in t direction
!! - Array dimensions of Fourier basis functions and derivatives
!!   (HZ, HZ_p):
!!   - Toroidal mode index
!!   - Toroidal plane index
!! - Array dimensions of one-dimensional basis functions
!!   (H1, H1_s, H1_ss)
!!   - Vertex
!!   - Basis function
!!   - Gaussian point
!!
!! @see Gauss
module basis_at_gaussian
  
  use mod_parameters
  use gauss
  
  implicit none
  
  save
  
  real*8 :: H(n_vertex_max,n_degrees,n_gauss,n_gauss)    !< Basis functions in poloidal plane
  real*8 :: H_s(n_vertex_max,n_degrees,n_gauss,n_gauss)  !< Basis functions in poloidal plane; derivative with respect to first coordinate
  real*8 :: H_t(n_vertex_max,n_degrees,n_gauss,n_gauss)  !< Basis functions in poloidal plane; derivative with respect to second coordinate
  real*8 :: H_st(n_vertex_max,n_degrees,n_gauss,n_gauss) !< Basis functions in poloidal plane; cross-derivative 
  real*8 :: H_ss(n_vertex_max,n_degrees,n_gauss,n_gauss) !< Basis functions in poloidal plane; second derivative with respect to first coordinate
  real*8 :: H_tt(n_vertex_max,n_degrees,n_gauss,n_gauss) !< Basis functions in poloidal plane; second derivative with respect to second coordinate

  real*8 :: HZ(n_tor,n_plane)          !< Basis functions in toroidal direction
  real*8 :: HZ_p(n_tor,n_plane)        !< Derivative of basis functions in toroidal direction
  real*8 :: HZ_pp(n_tor,n_plane)       !< Second derivative of basis functions in toroidal direction

  real*8 :: H1(2,n_degrees_1d,n_gauss)    !< One dimensional basis functions
  real*8 :: H1_s(2,n_degrees_1d,n_gauss)  !< First derivative of one dimensional basis functions
  real*8 :: H1_ss(2,n_degrees_1d,n_gauss) !< Second derivative of one dimensional basis functions
  
contains

!> calculates the basis functions at the Gaussian points
subroutine initialise_basis()
  use constants, only: PI
  use gauss
  use phys_module, only: n_tor, n_plane, n_period, mode
  use mod_basisfunctions

  implicit none

  ! --- local variables
  integer :: i,k,l
  real*8  :: s,t,phi

  do k=1,n_gauss
   s = xgauss(k)
   call basisfunctions1(s,H1(1:2,1:n_degrees_1d,k), H1_s(1:2,1:n_degrees_1d,k), H1_ss(1:2,1:n_degrees_1d,k)) ! the one-D basis functions
   do l=1,n_gauss
     t = xgauss(l)
     call basisfunctions(s,t,H(1:4,1:n_degrees,k,l),   H_s(1:4,1:n_degrees,k,l), H_t(1:4,1:n_degrees,k,l), &
                             H_st(1:4,1:n_degrees,k,l),H_ss(1:4,1:n_degrees,k,l),H_tt(1:4,1:n_degrees,k,l) )
   enddo
  enddo


  do k=1,n_plane

    phi = 2.d0*PI*float(k-1)/float(n_plane) / float(n_period)

    HZ(1,k)   = 1.d0
    HZ_p(1,k) = 0.d0

    do i=1,(n_tor-1)/2

      HZ(2*i,k)      =                           cos(mode(2*i)  *phi)
      HZ_p(2*i,k)    = - float(mode(2*i))      * sin(mode(2*i)  *phi)
      HZ_pp(2*i,k)   = - float(mode(2*i))**2   * cos(mode(2*i)  *phi)
      HZ(2*i+1,k)    =                           sin(mode(2*i+1)*phi)
      HZ_p(2*i+1,k)  = + float(mode(2*i+1))    * cos(mode(2*i+1)*phi)
      HZ_pp(2*i+1,k) = - float(mode(2*i+1))**2 * sin(mode(2*i+1)*phi)

    enddo

  enddo
end subroutine initialise_basis
end module basis_at_gaussian
