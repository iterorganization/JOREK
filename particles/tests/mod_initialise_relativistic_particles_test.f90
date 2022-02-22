!> test the initialisation routines for relativistic particles
module mod_initialise_relativistic_particles_test
use fruit
use constants,             only: PI,TWOPI
use mod_particle_sim,      only: particle_group
use mod_particle_types,    only: particle_kinetic_relativistic_id
use mod_particle_types,    only: particle_gc_relativistic_id
use mod_particle_types,    only: particle_kinetic_relativistic
use mod_particle_types,    only: particle_gc_relativistic
use mod_fields_analytical, only: fields_analytical
implicit none

private
public :: run_fruit_initialise_relativistic_particles

!> Variables and datatypes ------------------------------------
integer,parameter  :: n_groups=2
integer,parameter  :: n_samples=343
integer,dimension(n_groups),parameter    :: n_particles=(/123,234/)
integer,dimension(n_groups),parameter    :: p_types_sol=(/&
                      particle_kinetic_relativistic_id,&
                      particle_gc_relativistic_id/)
real*8,parameter                         :: tol_real8=5.d-13
real*8,parameter                         :: time_sol=0.d0
real*8,dimension(2),parameter            :: psi_big_sol=(/-2.d1,3.d2/)
real*8,dimension(2),parameter            :: psi_small_sol=(/0.d0,1.d0/)
real*8,dimension(2),parameter            :: psi_axisbnd=(/-1.d0,2.5d0/)
real*8,dimension(2),parameter            :: R_big_sol=(/1.d-1,3.d2/)
real*8,dimension(2),parameter            :: R_small_sol=(/1.d0,2.5d0/)
real*8,dimension(2),parameter            :: Z_big_sol=(/-3.d1,3.d1/)
real*8,dimension(2),parameter            :: Z_small_sol=(/-1.d0,2.5d0/)
real*8,dimension(2),parameter            :: theta_big_sol=(/-3.d1,3.d1/)
real*8,dimension(2),parameter            :: theta_small_sol=(/PI/6.d0,TWOPI/3.d0/)
real*8,dimension(2),parameter            :: phi_big_sol=(/-3.d1,3.d1/)
real*8,dimension(2),parameter            :: phi_small_sol=(/PI/6.d0,TWOPI/3.d0/)
real*8,dimension(2),parameter            :: psi_minmax=(/-5.d0,1.d1/)
real*8,dimension(2),parameter            :: R_minmax=(/5.d-1,1.d1/)
real*8,dimension(2),parameter            :: Z_minmax=(/-1.d1,1.d1/)
real*8,dimension(2),parameter            :: theta_minmax=(/0.d0,TWOPI/)
real*8,dimension(2),parameter            :: phi_minmax=(/0.d0,TWOPI/)
type(particle_group),dimension(n_groups) :: groups_sol
type(fields_analytical)                  :: fields_sol
!> Interfaces--------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> fruit test basket: set-up, run and tear-down tests
subroutine run_fruit_initialise_relativistic_particles()
  implicit none
  write(*,'(/A)') "  ... setting-up: initialise relativistic particles tests"
  call setup
  write(*,'(/A)') "  ... running: initialise relativistic particles tests"
  call test_dummy
  call test_sampling_cartesian_p_kinetic_relativistic
  call test_sampling_cartesian_gc_kinetic_relativistic
  call test_sampling_uniform_ppitchgyro_kinetic_relativistic
  call test_sampling_uniform_ppitchgyro_gc_relativistic
  call test_particle_base_init_to_zero
  call test_sampling_uniform_charge
  call test_check_psi_interval
  call test_check_RZPhi_interval
  write(*,'(/A)') "  ... tearing-down: initialise relativistic particles tests"
  call teardown
end subroutine run_fruit_initialise_relativistic_particles

!> Set-up and tear-down ---------------------------------------
!> set-up relativistic particle initialisation test features
subroutine setup()
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_particle_common_test_tools, only: fill_mass_RE
  use mod_particle_common_test_tools, only: rng_seed_interval
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  use mod_particle_common_test_tools, only: RZ0_lowbnd,RZ0_uppbnd
  use mod_particle_common_test_tools, only: BE0_lowbnd,BE0_uppbnd
  implicit none
  !> variables:
  integer :: ifail
  integer,dimension(0) :: int_param
  real*8,dimension(4)  :: real_param 
  !> allocate particle lists
  call allocate_one_particle_list_type(n_groups,n_particles,&
  p_types_sol,groups_sol,ifail)
  !> initialise the group masses as runaway
  call fill_mass_RE(n_groups,groups_sol)
  !> initialise particle fields
  call gnu_rng_interval(2,RZ0_lowbnd,RZ0_uppbnd,real_param(1:2))
  call gnu_rng_interval(2,BE0_lowbnd,BE0_uppbnd,real_param(3:4))
  call fields_sol%init_fields(0,4,int_param,real_param)
end subroutine setup

!> tear-down test feature
subroutine teardown()
  implicit none
  !> variables
  integer :: ii
  !> clean particle group
  do ii=1,n_groups; groups_sol(ii)%mass=0.d0; enddo
  !> cleanup fields
  call fields_sol%deallocate_fields()
end subroutine teardown

!> Tests ------------------------------------------------------
!> test initialisation particle momentum between limits
subroutine test_sampling_cartesian_p_kinetic_relativistic()
  use mod_gnu_rng,                           only: gnu_rng_interval
  use mod_particle_common_test_tools,        only: vp3d_lowbnd,vp3d_uppbnd
  use mod_initialise_relativistic_particles, only: sampling_cartesian_p_kinetic_relativistic
  implicit none
  !> variables:
  integer :: ii,jj
  real*8,dimension(3)  :: rand
  real*8,dimension(3,2) :: pxpypz_int
  logical,dimension(:),allocatable :: success
  !> initialisation
  pxpypz_int(:,1) = vp3d_lowbnd; pxpypz_int(:,2) = vp3d_uppbnd;
  !> the the random kinetic particle initialisation 
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_kinetic_relativistic)
      allocate(success(n_particles(ii))); success = .false.;
      do jj=1,n_particles(ii)
        call random_number(rand)
        call sampling_cartesian_p_kinetic_relativistic(p_list(jj),rand,pxpypz_int)
        success(jj) = &
        ((p_list(jj)%p(1).ge.pxpypz_int(1,1)).and.(p_list(jj)%p(1).le.pxpypz_int(1,2))).and.&
        ((p_list(jj)%p(2).ge.pxpypz_int(2,1)).and.(p_list(jj)%p(2).le.pxpypz_int(2,2))).and.&
        ((p_list(jj)%p(3).ge.pxpypz_int(3,1)).and.(p_list(jj)%p(3).le.pxpypz_int(3,2)))
        p_list(jj)%p =0.d0
      enddo
      call assert_true(all(success),"Error sampling cart. p kinetic relat.: momenta not in bound!")
      deallocate(success)
    end select
  enddo
end subroutine test_sampling_cartesian_p_kinetic_relativistic

!> test initialisation gc momentum between limits
subroutine test_sampling_cartesian_gc_kinetic_relativistic()
  use mod_gnu_rng,                           only: gnu_rng_interval
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_particle_common_test_tools,        only: RZPhi_lowbnd,RZPhi_uppbnd
  use mod_particle_common_test_tools,        only: vp3d_lowbnd,vp3d_uppbnd
  use mod_initialise_relativistic_particles, only: sampling_cartesian_p_gc_relativistic
  implicit none
  !> variables:
  integer :: ii,jj
  real*8 :: psi,U,B_norm
  real*8,dimension(2)   :: p_gc_sol
  real*8,dimension(3)   :: rand,B,E,RZPhi,p_kin_sol
  real*8,dimension(3,2) :: pxpypz_int
  logical,dimension(:),allocatable :: success_ppar,success_mu
  !> initialisation
  pxpypz_int(:,1) = vp3d_lowbnd; pxpypz_int(:,2) = vp3d_uppbnd;
  call gnu_rng_interval(3,RZPhi_lowbnd,RZPhi_uppbnd,RZPhi)
  call fields_sol%calc_EBPsiU(time_sol,0,RZPhi(1:2),RZPhi(3),&
  E,B,psi,U); B_norm = norm2(B); B = B/B_norm; 
  B = vector_cylindrical_to_cartesian(RZPhi(3),B)
  !> test the sampling of gc from cartesian momentum
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_gc_relativistic)
      allocate(success_ppar(n_particles(ii))); allocate(success_mu(n_particles(ii)));
      !> loops on the particles
      do jj=1,n_particles(ii)
        !> initialising particles
        p_list(jj)%x=RZPhi; p_list(jj)%st=RZPhi(1:2); p_list(jj)%i_elm=0;
        !> compute momentum
        call random_number(rand)
        call sampling_cartesian_p_gc_relativistic(p_list(jj),fields_sol,&
        rand,time_sol,groups_sol(ii)%mass,pxpypz_int)
        !> compute new momentum
        p_kin_sol = pxpypz_int(:,1)+(pxpypz_int(:,2)-pxpypz_int(:,1))*rand
        p_gc_sol(1) = dot_product(p_kin_sol,B);
        p_kin_sol = p_kin_sol - p_gc_sol(1)*B
        p_gc_sol(2) = dot_product(p_kin_sol,p_kin_sol)/(2.d0*groups_sol(ii)%mass*B_norm)
        !> checks
        success_ppar(jj) = .true.; success_mu(jj) = .true.;
        if(p_gc_sol(1).ne.0.d0) then
          success_ppar(jj) = (abs((p_gc_sol(1)-p_list(jj)%p(1))/p_gc_sol(1))).le.tol_real8
        else
          if(p_list(jj)%p(1).ne.0.d0) success_ppar(1) = .false.
        endif
        if(p_gc_sol(2).ne.0.d0) then 
          success_mu(jj) = (abs((p_gc_sol(2)-p_list(jj)%p(2))/p_gc_sol(2))).le.tol_real8
        else
          if(p_list(jj)%p(2).ne.0.d0) success_mu(jj) = .false.
        endif
        !> cleaning particles
        p_list(jj)%x=0.d0; p_list(jj)%st=0.d0; p_list(jj)%p=0.d0;
      enddo
      call assert_true(all(success_ppar),&
      "Error sampling cart. gc kinetic relat.: parallel momentum not in bound!")
      call assert_true(all(success_mu),&
      "Error sampling cart. gc kinetic relat.: magnetic moment not in bound!")
      deallocate(success_ppar); deallocate(success_mu);
    end select
  enddo
end subroutine test_sampling_cartesian_gc_kinetic_relativistic

!> test initialisation relativistic kinetic particle momentum from
!> momentum intensity, pitch and gyro angles
subroutine test_sampling_uniform_ppitchgyro_kinetic_relativistic()
  use constants,                             only: PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT,EL_CHG
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_gnu_rng,                           only: gnu_rng_interval
  use mod_pusher_tools,                      only: get_orthonormals
  use mod_particle_common_test_tools,        only: RZPhi_lowbnd,RZPhi_uppbnd
  use mod_particle_common_test_tools,        only: EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
  use mod_initialise_relativistic_particles, only: sampling_uniform_ppitchgyro_kinetic_relativistic
  implicit none
  !> variabes:
  integer :: ii,jj
  real*8 :: psi,U,B_norm,p_norm_test,theta_test,gyro_test
  real*8,dimension(2) :: p_int,costheta_int,gyro_int
  real*8,dimension(3) :: rand,RZPhi,B,E,e2,e3
  real*8,dimension(3,2) :: EThetaChi
  logical,dimension(:),allocatable :: success
  !> initialisations
  costheta_int= (/cos(EThetaChi_RE_lowbnd(2)),cos(EThetaChi_RE_uppbnd(2))/); 
  gyro_int = (/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/);
  call gnu_rng_interval(2,RZPhi_lowbnd,RZPhi_uppbnd,RZPhi)
  call fields_sol%calc_EBPsiU(time_sol,0,RZPhi(1:2),RZPhi(3),&
  E,B,psi,U); B_norm = norm2(B); B = B/B_norm; 
  B = vector_cylindrical_to_cartesian(RZPhi(3),B);
  !> compute orthogonal coordinate system
  call get_orthonormals(B,e2,e3)
  !> loop on the particles
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_kinetic_relativistic)
      p_int = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_int = (p_int+1.d0)*(p_int+1.0); p_int = sqrt(p_int-1.d0);
      p_int = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_int;
      allocate(success(n_particles(ii)))
      !> loop on the particles
      do jj=1,n_particles(ii)
        !> initialising particles
        call random_number(rand)
        p_list(jj)%x=RZPhi; p_list(jj)%st=RZPhi(1:2); p_list(jj)%i_elm=0;
        call sampling_uniform_ppitchgyro_kinetic_relativistic(p_list(jj),&
        fields_sol,rand,time_sol,p_int**3.d0,costheta_int,gyro_int)
        !> check solotion
        p_norm_test = norm2(p_list(jj)%p);
        theta_test = acos(dot_product(p_list(jj)%p,B)/p_norm_test)
        gyro_test  = PI+atan2(dot_product(p_list(jj)%p,e3),dot_product(p_list(jj)%p,e2))
        success(jj) = ((p_norm_test.ge.p_int(1)).and.(p_norm_test.le.p_int(2))).and.&
                      ((theta_test.ge.(EThetaChi_RE_lowbnd(2))).and.&
                      (theta_test.le.EThetaChi_RE_uppbnd(2))).and.&
                      ((gyro_test.ge.EThetaChi_RE_lowbnd(3)).and.&
                      (gyro_test.le.EThetaChi_RE_uppbnd(3)))
        !> cleaning particles
        p_list(jj)%x=0.d0; p_list(jj)%st=0.d0; p_list(jj)%p=0.d0;
      enddo
      call assert_true(all(success),&
      "Error sampling particle uniform p, pitch, gyro kinetic relat.: momenta not in bound!")
      deallocate(success)
    end select
  enddo
end subroutine test_sampling_uniform_ppitchgyro_kinetic_relativistic

!> test the initialisation of relativistic kinetic gcs from relativistic kinetic
!> particles using spherical coordinates
subroutine test_sampling_uniform_ppitchgyro_gc_relativistic()
  use constants,                             only: PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT,EL_CHG
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_gnu_rng,                           only: gnu_rng_interval
  use mod_particle_common_test_tools,        only: RZPhi_lowbnd,RZPhi_uppbnd
  use mod_particle_common_test_tools,        only: EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
  use mod_initialise_relativistic_particles, only: sampling_uniform_ppitchgyro_gc_relativistic
  implicit none
  !> variabes:
  integer :: ii,jj
  real*8 :: psi,U,B_norm,p_norm_test,theta_test,theta_test_2
  real*8,dimension(2) :: p_int,costheta_int,gyro_int
  real*8,dimension(3) :: rand,RZPhi,B,E,e2,e3
  real*8,dimension(3,2) :: EThetaChi
  logical,dimension(:),allocatable :: success
  !> initialisation
  costheta_int= (/cos(EThetaChi_RE_lowbnd(2)),cos(EThetaChi_RE_uppbnd(2))/); 
  gyro_int = (/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/);
  call gnu_rng_interval(2,RZPhi_lowbnd,RZPhi_uppbnd,RZPhi)
  call fields_sol%calc_EBPsiU(time_sol,0,RZPhi(1:2),RZPhi(3),&
  E,B,psi,U); B_norm = norm2(B); B = B/B_norm; 
  B = vector_cylindrical_to_cartesian(RZPhi(3),B);
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_gc_relativistic)
      p_int = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_int = (p_int+1.d0)*(p_int+1.0); p_int = sqrt(p_int-1.d0);
      p_int = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_int;
      allocate(success(n_particles(ii)))
      do jj=1,n_particles(ii)
        !> initialising particles
        call random_number(rand)
        p_list(jj)%x=RZPhi; p_list(jj)%st=RZPhi(1:2); p_list(jj)%i_elm=0;
        call sampling_uniform_ppitchgyro_gc_relativistic(p_list(jj),fields_sol,&
        rand,groups_sol(ii)%mass,time_sol,p_int**3.d0,costheta_int,gyro_int)
        !> checks
        p_norm_test = sqrt(p_list(jj)%p(1)*p_list(jj)%p(1) + &
        p_list(jj)%p(2)*2.d0*B_norm*groups_sol(ii)%mass) 
        theta_test = acos(p_list(jj)%p(1)/p_norm_test)
        theta_test_2 = asin(sqrt(p_list(jj)%p(2)*2.d0*B_norm*groups_sol(ii)%mass)/p_norm_test)
        success(jj) = ((p_norm_test.ge.p_int(1)).and.(p_norm_test.le.p_int(2))).and.&
                      ((theta_test.ge.(EThetaChi_RE_lowbnd(2))).and.&
                      (theta_test.le.EThetaChi_RE_uppbnd(2))).and.&
                      ((theta_test_2.ge.EThetaChi_RE_lowbnd(2)).and.&
                      (theta_test_2.le.EThetaChi_RE_uppbnd(2)))
      enddo
      call assert_true(all(success),&
      "Error sampling gc uniform p, pitch, gyro kinetic relat.: momenta not in bound!")
      deallocate(success)
    end select
  enddo
end subroutine test_sampling_uniform_ppitchgyro_gc_relativistic

!> test initialisation of particles to zero
subroutine test_particle_base_init_to_zero()
  use mod_initialise_relativistic_particles, only: init_particle_base_to_zero
  implicit none
  !> variables:
  integer :: ii,jj
  integer,dimension(:),allocatable  :: i_elm_sol,i_life_sol,i_elm_test,i_life_test
  real*4,dimension(:),allocatable   :: t_birth_sol,t_birth_test
  real*8,dimension(:),allocatable   :: weight_sol,weight_test
  real*8,dimension(:,:),allocatable :: x_sol,st_sol,x_test,st_test
  !> test initialisation to zero
  do ii=1,n_groups
    allocate(i_elm_sol(n_particles(ii)));   allocate(i_elm_test(n_particles(ii)));
    allocate(i_life_sol(n_particles(ii)));  allocate(i_life_test(n_particles(ii)));
    allocate(t_birth_sol(n_particles(ii))); allocate(t_birth_test(n_particles(ii)));
    allocate(weight_sol(n_particles(ii)));  allocate(weight_test(n_particles(ii)));
    allocate(x_sol(3,n_particles(ii)));     allocate(x_test(3,n_particles(ii)));
    allocate(st_sol(2,n_particles(ii)));      allocate(st_test(2,n_particles(ii)));
    i_elm_sol=0; i_life_sol=0; t_birth_sol=0.d0; weight_sol=1.d0; x_sol=0.d0; st_sol=0.d0;
    !> initialise particle to zero
    do jj=1,n_particles(ii)
      call init_particle_base_to_zero(groups_sol(ii)%particles(jj))
      i_elm_test(jj)   = groups_sol(ii)%particles(jj)%i_elm
      i_life_test(jj)  = groups_sol(ii)%particles(jj)%i_life
      t_birth_test(jj) = groups_sol(ii)%particles(jj)%t_birth
      weight_test(jj)  = groups_sol(ii)%particles(jj)%weight
      x_test(:,jj)     = groups_sol(ii)%particles(jj)%x
      st_test(:,jj)    = groups_sol(ii)%particles(jj)%st
    enddo
    !> checks
    call assert_equals(i_elm_test,i_elm_sol,n_particles(ii),&
    "initialise particle base to zero: mismatch!")
    call assert_equals(i_life_test,i_life_sol,n_particles(ii),&
    "initialise particle base to zero: mismatch!")
    call assert_equals(t_birth_test,t_birth_sol,n_particles(ii),&
    "initialise particle base to zero: mismatch!")
    call assert_equals(weight_test,weight_sol,n_particles(ii),&
    "initialise particle base to zero: mismatch!")
    call assert_equals(x_test,x_sol,3,n_particles(ii),&
    "initialise particle base to zero: x mismatch!")
    call assert_equals(st_test,st_sol,2,n_particles(ii),&
    "initialise particle base to zero: st mismatch!")
    !> cleanups
    deallocate(i_elm_sol);   deallocate(i_elm_test);
    deallocate(i_life_sol);  deallocate(i_life_test);
    deallocate(t_birth_sol); deallocate(t_birth_test);
    deallocate(weight_sol);  deallocate(weight_test);
    deallocate(x_sol);       deallocate(x_test);
    deallocate(st_sol);      deallocate(st_test);
  enddo
end subroutine test_particle_base_init_to_zero

!> test uniform random sampling
subroutine test_sampling_uniform_charge()
  use mod_particle_common_test_tools,        only: q1_pos_interval
  use mod_particle_common_test_tools,        only: q1_neg_interval
  use mod_particle_common_test_tools,        only: q1_posneg_interval
  use mod_initialise_relativistic_particles, only: sampling_uniform_charge
  implicit none
  !> variables
  integer                      :: ii,jj
  integer*1                    :: q1_test
  integer*1,dimension(2,5)     :: q1_intervals
  logical,dimension(n_samples) :: success
  real*8                       :: rand
  !> initialisation
  q1_intervals(:,1)=(/1,1/); q1_intervals(:,2)=(/-1,-1/);
  q1_intervals(:,3)=q1_pos_interval; q1_intervals(:,4)=q1_neg_interval;
  q1_intervals(:,5) = q1_posneg_interval;
  !> test charge sampling
  do ii=1,size(q1_intervals,2)
    do jj=1,n_samples
      call random_number(rand)
      q1_test = sampling_uniform_charge(rand,q1_intervals(:,ii))
      success(jj) = (q1_test.ge.q1_intervals(1,ii)).and.(q1_test.le.q1_intervals(2,ii))
    enddo
    !> checks
    call assert_true(all(success),"Error sampling uniform charge: charges not in interval!")
  enddo
end subroutine test_sampling_uniform_charge

!> test min max and renormalization of psi bound
subroutine test_check_psi_interval()
  use mod_initialise_relativistic_particles, only: check_psithetaphi_interval
  implicit none
  !> variables:
  real*8,dimension(2) :: psi_test,theta_test,phi_test
  !> check complete inflow
  psi_test=psi_small_sol; theta_test=theta_small_sol; phi_test=phi_small_sol;
  call check_psithetaphi_interval(psi_test,theta_test,phi_test,psi_axisbnd(1),&
  psi_axisbnd(2),psi_minmax,theta_minmax,phi_minmax)
  call assert_equals(psi_test,psi_axisbnd,2,&
  "Error check psi-theta-phi interval: psi inflow test mismatch!")
  call assert_equals(theta_test,theta_small_sol,2,&
  "Error check psi-theta-phi interval: theta inflow test mismatch!")
  call assert_equals(phi_test,phi_small_sol,2,&
  "Error check psi-theta-phi interval: phi inflow test mismatch!")
  !> check complete overflow
  psi_test=psi_big_sol; theta_test=theta_big_sol; phi_test=phi_big_sol;
  call check_psithetaphi_interval(psi_test,theta_test,phi_test,&
  psi_axisbnd(1),psi_axisbnd(2),psi_minmax,theta_minmax,phi_minmax)
  call assert_equals(psi_test,psi_minmax,2,&
  "Error check psi-theta-phi interval: psi overflow test mismatch!")
  call assert_equals(theta_test,theta_minmax,2,&
  "Error check psi-theta-phi interval: theta overflow test mismatch!")
  call assert_equals(phi_test,phi_minmax,2,&
  "Error check psi-theta-phi interval: phi overflow test mismatch!")
  !> checke underflow min and inflow max
  psi_test=(/psi_big_sol(1),psi_small_sol(2)/)
  theta_test=(/theta_big_sol(1),theta_small_sol(2)/)
  phi_test=(/phi_big_sol(1),phi_small_sol(2)/)
  call check_psithetaphi_interval(psi_test,theta_test,phi_test,&
  psi_axisbnd(1),psi_axisbnd(2),psi_minmax,theta_minmax,phi_minmax)
  call assert_equals(psi_test,(/psi_minmax(1),psi_axisbnd(2)/),2,&
  "Error check psi-theta-phi interval: psi min underflow test mismatch!")
  call assert_equals(theta_test,(/theta_minmax(1),theta_small_sol(2)/),2,&
  "Error check psi-theta-phi interval: theta min underflow test mismatch!")
  call assert_equals(phi_test,(/phi_minmax(1),phi_small_sol(2)/),2,&
  "Error check psi-theta-phi interval: phi min underflow test mismatch!")
  !> check inflow min and overflow max
  psi_test=(/psi_small_sol(1),psi_big_sol(2)/)
  theta_test=(/theta_small_sol(1),theta_big_sol(2)/)
  phi_test=(/phi_small_sol(1),phi_big_sol(2)/)
  call check_psithetaphi_interval(psi_test,theta_test,phi_test,&
  psi_axisbnd(1),psi_axisbnd(2),psi_minmax,theta_minmax,phi_minmax)
  call assert_equals(psi_test,(/psi_axisbnd(1),psi_minmax(2)/),2,&
  "Error check psi-theta-phi interval: psi max overflow test mismatch!")
  call assert_equals(theta_test,(/theta_small_sol(1),theta_minmax(2)/),2,&
  "Error check psi-theta-phi interval: theta max overflow test mismatch!")
  call assert_equals(phi_test,(/phi_small_sol(1),phi_minmax(2)/),2,&
  "Error check psi-theta-phi interval: phi max overflow test mismatch!")
end subroutine test_check_psi_interval

!> test RZPhi min max bounding
subroutine test_check_RZPhi_interval()
  use mod_initialise_relativistic_particles, only: check_RZPhi_interval
  implicit none
  !> variables
  real*8,dimension(2) :: R_test,Z_test,phi_test
  !> check complete inflow
  R_test=R_small_sol;Z_test=Z_small_sol;phi_test=phi_small_sol;
  call check_RZPhi_interval(R_test,Z_test,phi_test,&
  R_minmax,Z_minmax,phi_minmax)
  call assert_equals(R_test,R_small_sol,2,&
  "Error check RZPhi interval: R inflow test mismatch!")
  call assert_equals(Z_test,Z_small_sol,2,&
  "Error check RZPhi interval: Z inflow test mismatch!")
  call assert_equals(phi_test,phi_small_sol,2,&
  "Error check RZPhi interval: phi inflow test mismatch!")
  !> check complete overflow
  R_test=R_big_sol;Z_test=Z_big_sol;phi_test=phi_big_sol;
  call check_RZPhi_interval(R_test,Z_test,phi_test,&
  R_minmax,Z_minmax,phi_minmax)
  call assert_equals(R_test,R_minmax,2,&
  "Error check RZPhi interval: R overflow test mismatch!")
  call assert_equals(Z_test,Z_minmax,2,&
  "Error check RZPhi interval: Z overflow test mismatch!")
  call assert_equals(phi_test,phi_minmax,2,&
  "Error check RZPhi interval: phi overflow test mismatch!")
  !> check min underflow and max inflow
  R_test=(/R_big_sol(1),R_small_sol(2)/)
  Z_test=(/Z_big_sol(1),Z_small_sol(2)/)
  phi_test=(/phi_big_sol(1),phi_small_sol(2)/)
  call check_RZPhi_interval(R_test,Z_test,phi_test,&
  R_minmax,Z_minmax,phi_minmax)
  call assert_equals(R_test,(/R_minmax(1),R_small_sol(2)/),2,&
  "Error check RZPhi interval: R underflow test mismatch!")
  call assert_equals(Z_test,(/Z_minmax(1),Z_small_sol(2)/),2,&
  "Error check RZPhi interval: Z underflow test mismatch!")
  call assert_equals(phi_test,(/phi_minmax(1),phi_small_sol(2)/),2,&
  "Error check RZPhi interval: phi underflow test mismatch!")
  !> check min inflow and max overflow
  R_test=(/R_small_sol(1),R_big_sol(2)/)
  Z_test=(/Z_small_sol(1),Z_big_sol(2)/)
  phi_test=(/phi_small_sol(1),phi_big_sol(2)/)
  call check_RZPhi_interval(R_test,Z_test,phi_test,&
  R_minmax,Z_minmax,phi_minmax)
  call assert_equals(R_test,(/R_small_sol(1),R_minmax(2)/),2,&
  "Error check RZPhi interval: R overflow test mismatch!")
  call assert_equals(Z_test,(/Z_small_sol(1),Z_minmax(2)/),2,&
  "Error check RZPhi interval: Z overflow test mismatch!")
  call assert_equals(phi_test,(/phi_small_sol(1),phi_minmax(2)/),2,&
  "Error check RZPhi interval: phi overflow test mismatch!")
end subroutine test_check_RZPhi_interval

!> dummy test procedure
subroutine test_dummy()
  use mod_initialise_relativistic_particles
  use mod_fields_analytical, only: fields_analytical
  implicit none
  logical :: success
  success=.true.
  call assert_true(success,"Dummy test")
end subroutine test_dummy

!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_initialise_relativistic_particles_test
