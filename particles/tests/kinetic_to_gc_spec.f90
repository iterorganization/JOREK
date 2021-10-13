!> This module contains some tests verifying the transformation from kinetic to GC coordinates and vice versa
module kinetic_to_gc_spec
use mod_boris
use mod_particle_types
use data_structure
use fruit
implicit none

type(type_node_list)    :: node_list
type(type_element_list) :: element_list
real*8, parameter :: B(3) = [0d0, 0d0, 6d0]
real*8, parameter :: mass = 100 ! amu
real*8, parameter :: tol  = 1d-11
contains

subroutine setup_kinetic_to_gc_spec
  use mod_element_rtree, only: populate_element_rtree
  node_list%n_nodes = 0
  element_list%n_elements = 0

  call populate_element_rtree(node_list, element_list)
end subroutine setup_kinetic_to_gc_spec

!> convert a particle orbiting at 6 different points around a gc and see if we get the right answers
subroutine test_known_gc_kinetic
  use constants, only: PI, ATOMIC_MASS_UNIT, EL_CHG
  use mod_pusher_tools, only: get_orthonormals
  type(particle_kinetic) :: kinetic
  type(particle_gc)               :: gc

  real*8  :: chi, r, v_perp, v_par, e1(3), e2(3), x_ref(3), v_ref(3)
  integer :: i

  gc%x  = [1d0, 2d0, 3d0]
  gc%E  = 1d3 ! [eV]
  gc%mu = 1d2 ! [eV/T]
  gc%q  = 1
  gc%i_elm = 0

  ! gyroradius calculation
  v_perp = sqrt(2d0*gc%mu*EL_CHG*norm2(B)/(mass*ATOMIC_MASS_UNIT))
  v_par  = sqrt(2d0*(gc%E-gc%mu*norm2(B))*EL_CHG/(mass*ATOMIC_MASS_UNIT))
  r      = mass*ATOMIC_MASS_UNIT*v_perp/(gc%q*EL_CHG*(norm2(B)))

  ! check we did the calculation right
  call assert_equals(gc%E, (v_perp**2+v_par**2)*0.5d0*mass*ATOMIC_MASS_UNIT/EL_CHG,tol, "energy must be correct")

  call get_orthonormals(B, e1, e2)

  do i=1,12
    chi = PI*i/6d0
    kinetic = gc_to_kinetic(node_list, element_list, gc, chi, B, mass)

    x_ref = gc%x + cos(chi)*r*e2 - sin(chi)*r*e1
    v_ref = v_perp*cos(chi)*e1 + v_perp*sin(chi)*e2 + v_par*B/norm2(B)
    ! this uses the value of B being in the z-direction!
    call assert_equals(kinetic%x(1), x_ref(1), tol, "R positions must be same")
    call assert_equals(kinetic%x(2), x_ref(2), tol, "Z positions must be same")
    call assert_equals(kinetic%x(3), x_ref(3), tol, "phi positions must be same")
    call assert_equals(kinetic%v(1), v_ref(1), tol, "R velocity must be same")
    call assert_equals(kinetic%v(2), v_ref(2), tol, "Z velocity must be same")
    call assert_equals(kinetic%v(3), v_ref(3), tol, "phi velocity must be same")
  end do

end subroutine test_known_gc_kinetic

!> Test a few positions around a known guiding center and see if it is found okay
subroutine test_known_kinetic_gc
  use constants, only: PI, ATOMIC_MASS_UNIT, EL_CHG
  use mod_pusher_tools, only: get_orthonormals
  type(particle_kinetic) :: kinetic
  type(particle_gc)               :: gc

  real*8  :: chi, r, v_perp, v_par, e1(3), e2(3)
  integer :: i

  ! should do this in a prefilter
  node_list%n_nodes = 0
  element_list%n_elements = 0

  v_perp = 100.d0
  v_par  = 200.d0

  call get_orthonormals(B, e1, e2)
  do i=1,12
    chi = PI*i/6d0

    kinetic%q = 1
    kinetic%v = cos(chi)*v_perp*e1 + sin(chi)*v_perp*e2 + v_par*B/norm2(B) ! [m/s]
    r         = mass*ATOMIC_MASS_UNIT*norm2(kinetic%v(1:2))/(norm2(B)*kinetic%q*EL_CHG) ! m v_perp / qB
    kinetic%x = [1.d0, 2.d0, 3.d0] + r * (cos(chi)*e2-sin(chi)*e1)
    kinetic%i_elm = 0

    gc = kinetic_to_gc(node_list, element_list, kinetic, B, mass)
    ! this uses the value of B being in the z-direction!
    call assert_equals(1d0, gc%x(1), tol, "R positions must be same")
    call assert_equals(2d0, gc%x(2), tol, "Z positions must be same")
    call assert_equals(3d0, gc%x(3), tol, "phi positions must be same")
    call assert_equals(gc%E,  0.5d0*mass*ATOMIC_MASS_UNIT*dot_product(kinetic%v,kinetic%v)/EL_CHG, tol, "Energy must be right")
    call assert_equals(gc%mu, 0.5d0*mass*ATOMIC_MASS_UNIT*dot_product(kinetic%v(1:2),kinetic%v(1:2))/norm2(B)/EL_CHG, tol, "mu must be right")
  end do

end subroutine test_known_kinetic_gc

!> Test the performance by transforming back and forth
subroutine test_kinetic_gc_kinetic
  use mod_pusher_tools, only: get_orthonormals
  use constants, only: PI
  type(particle_kinetic) :: kinetic1, kinetic2
  type(particle_gc)               :: gc
  real*8 :: chi, e1(3), e2(3)

  kinetic1%x  = [1d0, 2d0, 3d0]
  kinetic1%v  = [1d3, 1d2, 2d2]
  kinetic1%q  = 1
  kinetic1%i_elm = 0
  ! angle of velocity vector with r-axis + pi to get angle of gyration with r-axis
  call get_orthonormals(B, e1, e2)
  chi = atan2(dot_product(kinetic1%v,e2),dot_product(kinetic1%v,e1))

  gc = kinetic_to_gc(node_list, element_list, kinetic1, B, mass)
  kinetic2 = gc_to_kinetic(node_list, element_list, gc, chi, B, mass)

  call assert_equals(kinetic1%x(1), kinetic2%x(1), tol, "R positions must be same")
  call assert_equals(kinetic1%x(2), kinetic2%x(2), tol, "Z positions must be same")
  call assert_equals(kinetic1%x(3), kinetic2%x(3), tol, "phi positions must be same")
  call assert_equals(kinetic1%v(1), kinetic2%v(1), tol, "R velocity must be same")
  call assert_equals(kinetic1%v(2), kinetic2%v(2), tol, "Z velocity must be same")
  call assert_equals(kinetic1%v(3), kinetic2%v(3), tol, "phi velocity must be same")
  call assert_equals(int(kinetic1%q,4), int(kinetic2%q,4), "q must be same")
  ! Do not check st, i_elm
end subroutine test_kinetic_gc_kinetic

!> Test the performance by transforming back and forth
subroutine test_gc_kinetic_gc
  type(particle_kinetic) :: kinetic
  type(particle_gc)               :: gc1, gc2
  real*8, parameter :: chi  = 0.d0

  gc1%x  = [1d0, 1d0, 2d0]
  gc1%E  = 1d3
  gc1%mu = 1d2
  gc1%q  = 1
  gc1%i_elm = 0

  kinetic = gc_to_kinetic(node_list, element_list, gc1, chi, B, mass)
  gc2 = kinetic_to_gc(node_list, element_list, kinetic, B, mass)

  call assert_equals(gc1%x(1), gc2%x(1), tol, "R positions must be same")
  call assert_equals(gc1%x(2), gc2%x(2), tol, "Z positions must be same")
  call assert_equals(gc1%x(3), gc2%x(3), tol, "phi positions must be same")
  call assert_equals(gc1%E, gc2%E, tol, "Energies must be same")
  call assert_equals(gc1%mu, gc2%mu, tol, "Mu must be same")
  call assert_equals(int(gc1%q,4), int(gc2%q,4), "q must be same")
  ! Do not check st, i_elm
end subroutine test_gc_kinetic_gc

!> Test the performance by transforming back and forth
subroutine test_gc_kinetic_gc_negative_mu
  type(particle_kinetic) :: kinetic
  type(particle_gc)               :: gc1, gc2
  real*8, parameter :: chi  = 0.d0

  gc1%x  = [1d0, 1d0, 2d0]
  gc1%E  = 1d3
  gc1%mu = -1d2
  gc1%q  = 1
  gc1%i_elm = 0

  kinetic = gc_to_kinetic(node_list, element_list, gc1, chi, B, mass)
  gc2 = kinetic_to_gc(node_list, element_list, kinetic, B, mass)

  call assert_equals(gc1%x(1), gc2%x(1), tol, "R positions must be same")
  call assert_equals(gc1%x(2), gc2%x(2), tol, "Z positions must be same")
  call assert_equals(gc1%x(3), gc2%x(3), tol, "phi positions must be same")
  call assert_equals(gc1%E, gc2%E, tol, "Energies must be same")
  call assert_equals(gc1%mu, gc2%mu, tol, "Mu must be same")
  call assert_equals(int(gc1%q,4), int(gc2%q,4), "q must be same")
  ! Do not check st, i_elm
end subroutine test_gc_kinetic_gc_negative_mu

!> Test function to generate orthonormal vectors
subroutine test_get_orthonormals
  use mod_pusher_tools, only: get_orthonormals

  real*8, dimension(3) :: in, e1, e2

  in = [1d0, 0d0, 0d0]
  call get_orthonormals(in, e1, e2)
  call do_test_orthonormality(in, e1, e2)

  in = [2d0, 0d0, 0d0]
  call get_orthonormals(in, e1, e2)
  call do_test_orthonormality(in, e1, e2)

  in = [-2d0, 0d0, 0d0]
  call get_orthonormals(in, e1, e2)
  call do_test_orthonormality(in, e1, e2)

  in = [0d0, 1d0, 0d0]
  call get_orthonormals(in, e1, e2)
  call do_test_orthonormality(in, e1, e2)

  in = [0d0, 1d0, 0d0]
  call get_orthonormals(in, e1, e2)
  call do_test_orthonormality(in, e1, e2)

  in = [0d0, -2d0, 0d0]
  call get_orthonormals(in, e1, e2)
  call do_test_orthonormality(in, e1, e2)

  in = [-2d0, 5d0, -1d0]
  call get_orthonormals(in, e1, e2)
  call do_test_orthonormality(in, e1, e2)

  in = [-2d0, -5d0, -1d0]
  call get_orthonormals(in, e1, e2)
  call do_test_orthonormality(in, e1, e2)
end subroutine test_get_orthonormals

subroutine do_test_orthonormality(in, e1, e2)
  real*8, dimension(3), intent(in) :: in, e1, e2

  call assert_equals(norm2(e1), 1.d0, tol, "Norm of 1 must be 1")
  call assert_equals(norm2(e2), 1.d0, tol, "Norm of 2 must be 1")
  call assert_equals(dot_product(in,e1), 0.d0, tol, "1 must be orthogonal to in")
  call assert_equals(dot_product(in,e2), 0.d0, tol, "2 must be orthogonal to in")
  call assert_equals(dot_product(e1,e2), 0.d0, tol, "1 must be orthogonal to 2")
end subroutine do_test_orthonormality
end module kinetic_to_gc_spec
