module mod_kinetic_to_gc_spec_test
use fruit
use mod_boris
use mod_particle_types
use data_structure
implicit none
private
public :: run_fruit_kinetic_to_gc_spec

!> Variables --------------------------------------
type(type_node_list)    :: node_list
type(type_element_list) :: element_list
real*8, parameter :: B(3) = [0d0, 0d0, 6d0]
real*8, parameter :: mass = 100 ! amu
real*8, parameter :: tol  = 1d-11

contains
!> Fruit basket -----------------------------------
subroutine run_fruit_kinetic_to_gc_spec
  implicit none
  call init_node_list(node_list, n_nodes_max, node_list%n_dof, n_var)
  write(*,'(/A)') "  ... setting-up: kinetic to gc spec"
  call setup_kinetic_to_gc_spec
  write(*,'(/A)') "  ... running: kinetic to gc spec"
  call run_test_case(test_known_gc_kinetic,'test_known_gc_kinetic')
  call run_test_case(test_known_kinetic_gc,'test_known_kinetic_gc')
  call run_test_case(test_kinetic_gc_kinetic,'test_kinetic_gc_kinetic')
  call run_test_case(test_gc_kinetic_gc,'test_gc_kinetic_gc')
  call run_test_case(test_gc_kinetic_gc_negative_mu,'test_gc_kinetic_gc_negative_mu')
  call run_test_case(test_get_orthonormals,'test_get_orthonormals')
  call run_test_case(test_cartesian_conversions,'test_cartesian_conversions')
  write(*,'(/A)') "  ... tearing-down: kinetic to gc spec"
  call dealloc_node_list(node_list)
end subroutine run_fruit_kinetic_to_gc_spec

!> Set-up and tear-down  --------------------------
subroutine setup_kinetic_to_gc_spec
  use mod_element_rtree, only: populate_element_rtree
  node_list%n_nodes = 0
  element_list%n_elements = 0

  call populate_element_rtree(node_list, element_list)
end subroutine setup_kinetic_to_gc_spec

!> Tests ------------------------------------------

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

!> Tools ------------------------------------------
subroutine do_test_orthonormality(in, e1, e2)
  real*8, dimension(3), intent(in) :: in, e1, e2
  call assert_equals(norm2(e1), 1.d0, tol, "Norm of 1 must be 1")
  call assert_equals(norm2(e2), 1.d0, tol, "Norm of 2 must be 1")
  call assert_equals(dot_product(in,e1), 0.d0, tol, "1 must be orthogonal to in")
  call assert_equals(dot_product(in,e2), 0.d0, tol, "2 must be orthogonal to in")
  call assert_equals(dot_product(e1,e2), 0.d0, tol, "1 must be orthogonal to 2")
end subroutine do_test_orthonormality
!> ------------------------------------------------
! A uniform CARTESIAN field must be expressed at each particle's local phi.
! Independent Cartesian formulas below deliberately do not use the geometry helper.
subroutine test_cartesian_conversions
  use constants, only: ATOMIC_MASS_UNIT, EL_CHG, TWOPI, SPEED_OF_LIGHT
  use mod_pusher_tools, only: get_orthonormals
  use mod_math_operators, only: cross_product
  use mod_gc_variational, only: convert_gc_vpar_to_kinetic, convert_leapfrog_to_gc_vpar
  use mod_gc_relativistic, only: relativistic_gc_to_relativistic_kinetic
  use mod_kinetic_relativistic, only: gc_to_relativistic_kinetic, &
      relativistic_kinetic_to_gc, relativistic_kinetic_to_relativistic_gc
  type(particle_gc) :: gc, back
  type(particle_gc_vpar) :: gv, gv_back
  type(particle_gc_relativistic) :: gr, gr_back
  type(particle_kinetic) :: kin
  type(particle_kinetic_leapfrog) :: orbit(4)
  type(particle_kinetic_relativistic) :: kr
  real*8 :: field(3), bh(3), e1(3), e2(3), v(3), vc(3), bc(3), rho(3), bp(3), pc(3)
  real*8 :: radii(3), angles(3), dirs(3,3), bn, phase, u, energy, moment, test_mass
  integer :: i,j,k,q,fail,nseed
  integer, allocatable :: saved(:), seed(:), expected(:), actual(:)

  radii=[.23d0,1.7d0,4.6d0]; angles=[-.2d0,TWOPI-1.d-4,3.d0*TWOPI+.2d0]
  dirs(:,1)=[0.d0,0.d0,1.d0];dirs(:,2)=[0.d0,1.d0,0.d0];dirs(:,3)=[.3d0,.8d0,.5d0]
  test_mass=2.01410178d0
  call random_seed(size=nseed);allocate(saved(nseed),seed(nseed),expected(nseed),actual(nseed))
  call random_seed(get=saved)
  do i=1,3
  do j=1,3
  do k=1,3
  do q=-1,1,2
    field=2.3d0*dirs(:,i)/norm2(dirs(:,i));bn=norm2(field);bh=field/bn
    gv%x=[radii(j),.13d0,angles(k)];gv%st=0.d0;gv%i_elm=0;gv%q=int(q,1);gv%weight=1.d0
    gv%vpar=q*1.7d5;gv%mu=2.1d10;gv%B_norm=bn
    ! Determine the existing randomized phase without changing its production API.
    seed=9173+13*i+7*j+3*k
    call random_seed(put=seed);call random_number(u);phase=TWOPI*u
    call random_seed(get=expected);call random_seed(put=seed)
    call convert_gc_vpar_to_kinetic(node_list,element_list,gv,field,test_mass,4,orbit,fail)
    call random_seed(get=actual)
    call assert_true(all(actual==expected),'positive-N conversion consumes exactly one RNG draw')
    call get_orthonormals(bh,e1,e2)
    v=gv%vpar*bh+sqrt(2.d0*gv%mu*bn)*(cos(phase)*e1+sin(phase)*e2)
    vc=cart_vector(gv%x(3),v);bc=cart_vector(gv%x(3),bh)
    rho=test_mass*ATOMIC_MASS_UNIT/(q*EL_CHG*bn)*cart_cross(bc,vc)
    call assert_true(norm2(cart_vector(gv%x(3),-cross_product(v,bh))-cart_cross(bc,vc)) &
        <tol*norm2(v),'existing cross-product sign agrees with Cartesian determinant')
    call assert_true(norm2(cart_position(orbit(1)%x)-cart_position(gv%x)-rho)<tol, &
        'vpar conversion: Cartesian displacement, including toroidal length')
    call assert_true(norm2(cart_vector(orbit(1)%x(3),orbit(1)%v)-vc)<tol*norm2(v), &
        'vpar conversion: physical Cartesian velocity is unchanged')
    bp=local_vector(orbit(1)%x(3),bn*bc)
    call convert_leapfrog_to_gc_vpar(node_list,element_list,orbit(1),bp,test_mass,gv_back)
    call assert_true(norm2(cart_position(gv_back%x)-cart_position(gv%x))<tol,'vpar inverse position')
    call assert_equals(gv%x(3),gv_back%x(3),tol,'continuous phi round trip')
    call assert_equals(1.d0,gv_back%vpar/gv%vpar,tol,'vpar inverse parallel velocity')
    call assert_equals(1.d0,gv_back%mu/gv%mu,tol,'vpar inverse magnetic moment')
    call assert_equals(1.d0,sum(orbit(1)%v**2)/(gv%vpar**2+2.d0*gv%mu*bn),tol,'kinetic energy')

    energy=.5d0*test_mass*ATOMIC_MASS_UNIT*sum(v**2)/EL_CHG
    moment=sign(test_mass*ATOMIC_MASS_UNIT*gv%mu/EL_CHG,gv%vpar)
    gc=gv;gc%q=gv%q;gc%E=energy;gc%mu=moment
    kin=gc_to_kinetic(node_list,element_list,gc,phase,field,test_mass)
    call assert_true(norm2(cart_position(kin%x)-cart_position(gc%x)-rho)<tol,'E/mu forward position')
    call assert_true(norm2(cart_vector(kin%x(3),kin%v)-vc)<tol*norm2(v),'E/mu velocity basis')
    bp=local_vector(kin%x(3),bn*bc)
    back=kinetic_to_gc(node_list,element_list,kin,bp,test_mass)
    call assert_true(norm2(cart_position(back%x)-cart_position(gc%x))<tol,'E/mu inverse position')
    call assert_equals(1.d0,back%E/energy,tol,'E/mu inverse energy')
    call assert_equals(1.d0,back%mu/moment,tol,'E/mu inverse signed magnetic moment')

    ! Larger momenta avoid ill-conditioning from subtracting total rest energy
    ! in the legacy relativistic E/mu interface; no test tolerance is relaxed.
    gr=gv;gr%q=gv%q;gr%p=[100.d0*test_mass*gv%vpar,1.d4*test_mass*gv%mu]
    pc=100.d0*test_mass*vc
    kr=relativistic_gc_to_relativistic_kinetic(node_list,element_list,gr,test_mass,field,phase)
    call assert_true(norm2(cart_position(kr%x)-cart_position(gv%x)-100.d0*rho)<tol,'relativistic position')
    call assert_true(norm2(kr%p-pc)<tol*norm2(pc),'relativistic Cartesian momentum uses GC angle')
    bp=local_vector(kr%x(3),bn*bc)
    gr_back=relativistic_kinetic_to_relativistic_gc(node_list,element_list,kr,test_mass,bp)
    call assert_true(norm2(cart_position(gr_back%x)-cart_position(gr%x))<tol,'relativistic inverse position')
    call assert_true(all(abs(gr_back%p/gr%p-1.d0)<tol),'relativistic inverse momenta')
    gc%E=ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*sqrt((test_mass*SPEED_OF_LIGHT)**2+sum(pc**2))/EL_CHG
    gc%mu=moment*1.d4
    kr=gc_to_relativistic_kinetic(node_list,element_list,gc,0.d0,test_mass,phase,field)
    call assert_true(norm2(kr%p-pc)<tol*norm2(pc),'legacy relativistic Cartesian momentum')
    call assert_true(norm2(cart_position(kr%x)-cart_position(gc%x)-100.d0*rho)<tol, &
        'legacy relativistic forward position')
    bp=local_vector(kr%x(3),bn*bc)
    back=relativistic_kinetic_to_gc(node_list,element_list,kr,test_mass,bp)
    call assert_true(norm2(cart_position(back%x)-cart_position(gc%x))<tol,'legacy relativistic inverse position')
    call assert_equals(1.d0,back%E/gc%E,tol,'relativistic total energy')
    call assert_equals(1.d0,back%mu/gc%mu,tol,'relativistic signed magnetic moment')
  enddo
  enddo
  enddo
  enddo
  call random_seed(get=expected)
  call convert_gc_vpar_to_kinetic(node_list,element_list,gv,field,test_mass,0,orbit,fail)
  call random_seed(get=actual)
  call assert_true(all(actual==expected),'zero-N conversion consumes no RNG')
  call random_seed(put=saved)
end subroutine test_cartesian_conversions

pure function cart_position(x) result(y)
  real*8,intent(in)::x(3)
  real*8::y(3)
  y=[x(1)*cos(x(3)),-x(1)*sin(x(3)),x(2)]
end function
pure function cart_vector(phi,v) result(w)
  real*8,intent(in)::phi,v(3)
  real*8::w(3)
  w=[cos(phi)*v(1)-sin(phi)*v(3),-sin(phi)*v(1)-cos(phi)*v(3),v(2)]
end function
pure function local_vector(phi,v) result(w)
  real*8,intent(in)::phi,v(3)
  real*8::w(3)
  w=[cos(phi)*v(1)-sin(phi)*v(2),v(3),-sin(phi)*v(1)-cos(phi)*v(2)]
end function
pure function cart_cross(a,b) result(c)
  real*8,intent(in)::a(3),b(3)
  real*8::c(3)
  c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
end function

end module mod_kinetic_to_gc_spec_test
