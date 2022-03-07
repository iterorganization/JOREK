!> mod_fields_minmax_test implements variables and procedures
!> used for testing methods use for finding jorek fields
!> minima and maxima
module mod_fields_minmax_mpi_test
use fruit
use mod_rng,           only: type_rng
use mod_fields_minmax, only: fun_interp_PRZ
implicit none
private
public :: run_fruit_fields_minmax_mpi

!> Variables and datatypes ----------------------------------
character(len=31),parameter :: filename="test_jorek_3d_fields_restart.h5"
integer,parameter           :: rst_format=0
integer,parameter           :: n_coords=3
integer,parameter           :: n_s=37
integer,parameter           :: n_t=25
integer,parameter           :: n_phi=100
integer,parameter           :: field_id_sol=1 !< poloidal flux
real*8,parameter            :: tol_real8=1.d-15
type(fun_interp_PRZ)        :: interp_PRZ_object
class(type_rng),dimension(:),allocatable :: rngs

!> Interfaces -----------------------------------------------
contains 
!> Fruit basket ---------------------------------------------
!> wrapper for all test features set-up, run and tear-down procedures
subroutine run_fruit_fields_minmax_mpi(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  write(*,'(/A)') "  ... setting-up: fields minmax tests"
  call setup(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... running: fields minmax tests"
  call test_generation_equidistant_mesh_1d(rank,n_tasks,ifail)
  call test_generation_random_mesh_1d(rank,n_tasks,ifail)
  !call test_field_minmax_equidistant_mesh(rank,n_tasks,ifail)
  !call test_field_minmax_random_mesh(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... tearing-down: fields minmax tests"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_fields_minmax_mpi

!> Set-up and tear-down -------------------------------------
!> set-up the fields minmax features
subroutine setup(rank,n_tasks,ifail)
  use mod_pcg32_rng,      only: pcg32_rng
  use mod_rng,            only: setup_shared_rngs
  use mod_import_restart, only: import_hdf5_restart
  use basis_at_gaussian,  only: initialise_basis
  use data_structure,     only: init_threads
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> initialise the jorek simulation
  call init_threads
  call det_modes
  call initialise_and_broadcast_parameters(rank,"__NO_FILENAME__")
  call broadcast_phys(rank)
  call initialise_basis
  !> read jorek restart file
  call import_hdf5_restart(interp_PRZ_object%node_list,&
  interp_PRZ_object%element_list,trim(filename),rst_format,ifail)
  !> set-up the random number generator
  call setup_shared_rngs(n_coords,pcg32_rng(),rngs)
end subroutine setup

!> tear-down all test features
subroutine teardown(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> deallocate rng
  deallocate(rngs)
end subroutine teardown

!> Tests ----------------------------------------------------
!> test the generation of equidistant mesh in s,t,phi
subroutine test_generation_equidistant_mesh_1d(rank,n_tasks,ifail)
  use constants,         only: PI,TWOPI
  use mod_fields_minmax, only: generate_mesh_1d
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> variables
  integer :: ii
  real*8 :: d_phi_sol,d_st_sol
  real*8,dimension(n_phi)   :: eq_mesh_phi_test,eq_mesh_phi_sol
  real*8,dimension(n_s*n_t) :: eq_mesh_st_test,eq_mesh_st_sol
  !> initialisations
  d_st_sol = 1.d0/real(n_s*n_t-1,kind=8); d_phi_sol = PI/real(n_phi-1,kind=8);
  !> generate equidistant mesh
  call generate_mesh_1d(n_s*n_t,(/0.d0,1.d0/),eq_mesh_st_test)
  call generate_mesh_1d(n_phi,(/PI,TWOPI/),eq_mesh_phi_test)
  !> checks
  do ii=1,n_phi
    eq_mesh_phi_sol(ii) = PI+d_phi_sol*(ii-1)
  enddo
  do ii=1,n_s*n_t
    eq_mesh_st_sol(ii) = d_st_sol*(ii-1)
  enddo
  call assert_equals(eq_mesh_st_test,eq_mesh_st_sol,n_s*n_t,tol_real8,&
  "Error generate equidistant mesh 1d: mesh in (0,1) mismatch!")
  call assert_equals(eq_mesh_st_test(1),0.d0,tol_real8,&
  "Error generate equidistant mesh 1d: mesh in (0,1) lower bound mismatch!")
  call assert_equals(eq_mesh_st_test(n_s*n_t),1.d0,tol_real8,&
  "Error generate equidistant mesh 1d: mesh in (0,1) upper bound mismatch!")
  call assert_equals(eq_mesh_phi_test,eq_mesh_phi_sol,n_phi,tol_real8,&
  "Error generate equidistant mesh 1d: mesh in interval mismatch!")
  call assert_equals(eq_mesh_phi_test(1),PI,tol_real8,&
  "Error generate equidistant mesh 1d: mesh in interval lower bound mismatch!")
  call assert_equals(eq_mesh_phi_test(n_phi),TWOPI,tol_real8,&
  "Error generate equidistant mesh 1d: mesh in interval upper bound mismatch!")
end subroutine test_generation_equidistant_mesh_1d

!> test the generation of random mesh
subroutine test_generation_random_mesh_1d(rank,n_tasks,ifail)
  use constants,         only: PI,TWOPI
  use mod_fields_minmax, only: generate_mesh_1d
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> variables
  integer :: ii,jj
  real*8,dimension(n_phi)   :: rand_mesh_phi
  real*8,dimension(n_s*n_t) :: rand_mesh_st
  !> initialisation
  !> compute random mesh
  call generate_mesh_1d(n_s*n_t,rngs,(/0.d0,1.d0/),rand_mesh_st)
  call generate_mesh_1d(n_phi,rngs,(/PI,TWOPI/),rand_mesh_phi)
  !> checks
  call assert_true(all((rand_mesh_st.ge.0.d0).and.(rand_mesh_st.le.1.d0)),&
  "Error fields minmax random mesh 1d: nodes not in bound (0,1)!")
  call assert_true(all((rand_mesh_phi.ge.PI).and.(rand_mesh_phi.le.TWOPI)),&
  "Error fields minmax random mesh: nodes not in interval!")
end subroutine test_generation_random_mesh_1d

!> test the find of extrema on random mesh
subroutine test_field_minmax_random_mesh(rank,n_tasks,ifail)
  use constants,         only: TWOPI
  use mod_interp,        only: interp_PRZ
  use mod_fields_minmax, only: generate_mesh_1d
  use mod_fields_minmax, only: field_minmax
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> variables:
  integer             :: ii,jj,kk,pp
  real*8              :: R_test,Z_test
  real*8,dimension(1) :: field_test
  real*8,dimension(2) :: minmax_global
  real*8,dimension(interp_PRZ_object%element_list%n_elements) :: min_list,max_list
  real*8,dimension(n_s)   :: rand_mesh_s
  real*8,dimension(n_t)   :: rand_mesh_t
  real*8,dimension(n_phi) :: rand_mesh_phi
  logical,dimension(interp_PRZ_object%element_list%n_elements) :: fail_min_list
  logical,dimension(interp_PRZ_object%element_list%n_elements) :: fail_max_list
  logical,dimension(2) :: fail_global
  !> initialisation
  fail_min_list = .false.; fail_max_list = .false.; fail_global = .false.;
  call generate_mesh_1d(n_s,rngs,(/0.d0,1.d0/),rand_mesh_s)
  call generate_mesh_1d(n_t,rngs,(/0.d0,1.d0/),rand_mesh_t)
  call generate_mesh_1d(n_phi,rngs,(/0.d0,TWOPI/),rand_mesh_phi)
  !> find extrema
  call field_minmax(field_id_sol,n_s,n_t,n_phi,rand_mesh_s,rand_mesh_t,&
  rand_mesh_phi,interp_PRZ_object,min_list,max_list,minmax_global)
  !> checks if the extrema are the largest and smallest values for each point in the element
  do ii=1,interp_PRZ_object%element_list%n_elements
    do jj=1,n_phi
      do kk=1,n_t
        do pp=1,n_s
          call interp_PRZ(interp_PRZ_object%node_list,interp_PRZ_object%element_list,&
          ii,(/field_id_sol/),1,rand_mesh_s(pp),rand_mesh_t(kk),rand_mesh_phi(jj),&
          field_test,R_test,Z_test)
          if(.not.fail_min_list(ii)) fail_min_list(ii) = field_test(1).lt.min_list(ii)
          if(.not.fail_max_list(ii)) fail_max_list(ii) = field_test(1).gt.max_list(ii)
          if(.not.fail_global(1))    fail_global(1)    = field_test(1).lt.minmax_global(1)
          if(.not.fail_global(2))    fail_global(2)    = field_test(1).gt.minmax_global(2)
        enddo
      enddo
    enddo
  enddo
  call assert_true(all(.not.fail_min_list),&
  "Error find fields minmax random mesh: minimum list not a minimum!")
  call assert_true(all(.not.fail_max_list),&
  "Error find fields minmax random mesh: maximum list not a maximum!")
  call assert_true(all(.not.fail_global),&
  "Error find fields minmax random mesh: not global minimum or maximum!")
end subroutine test_field_minmax_random_mesh

!> test the find of extrema on equidistant mesh
subroutine test_field_minmax_equidistant_mesh(rank,n_tasks,ifail)
  use constants,         only: PI,TWOPI
  use mod_interp,        only: interp_PRZ
  use mod_fields_minmax, only: generate_mesh_1d
  use mod_fields_minmax, only: field_minmax
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> variables:
  integer             :: ii,jj,kk,pp
  real*8              :: R_test,Z_test
  real*8,dimension(1) :: field_test
  real*8,dimension(2) :: minmax_global
  real*8,dimension(interp_PRZ_object%element_list%n_elements) :: min_list,max_list
  real*8,dimension(n_s)   :: eq_mesh_s
  real*8,dimension(n_t)   :: eq_mesh_t
  real*8,dimension(n_phi) :: eq_mesh_phi
  logical,dimension(interp_PRZ_object%element_list%n_elements) :: fail_min_list
  logical,dimension(interp_PRZ_object%element_list%n_elements) :: fail_max_list
  logical,dimension(2) :: fail_global
  !> initialisation
  fail_min_list = .false.; fail_max_list = .false.; fail_global = .false.;
  call generate_mesh_1d(n_s,(/0.d0,1.d0/),eq_mesh_s)
  call generate_mesh_1d(n_t,(/0.d0,1.d0/),eq_mesh_t)
  call generate_mesh_1d(n_phi,(/0.d0,TWOPI/),eq_mesh_phi)
  !> find extrema
  call field_minmax(field_id_sol,n_s,n_t,n_phi,eq_mesh_s,eq_mesh_t,&
  eq_mesh_phi,interp_PRZ_object,min_list,max_list,minmax_global)
  !> checks if the extrema are the largest and smallest values for each point in the element
  do ii=1,interp_PRZ_object%element_list%n_elements
    do jj=1,n_phi
      do kk=1,n_t
        do pp=1,n_s
          call interp_PRZ(interp_PRZ_object%node_list,interp_PRZ_object%element_list,&
          ii,(/field_id_sol/),1,eq_mesh_s(pp),eq_mesh_t(kk),eq_mesh_phi(jj),&
          field_test,R_test,Z_test)
          if(.not.fail_min_list(ii)) fail_min_list(ii) = field_test(1).lt.min_list(ii)
          if(.not.fail_max_list(ii)) fail_max_list(ii) = field_test(1).gt.max_list(ii)
          if(.not.fail_global(1))    fail_global(1)    = field_test(1).lt.minmax_global(1)
          if(.not.fail_global(2))    fail_global(2)    = field_test(1).gt.minmax_global(2)
        enddo
      enddo
    enddo
  enddo
  call assert_true(all(.not.fail_min_list),&
  "Error find fields minmax equidistant mesh: minimum list not a minimum!")
  call assert_true(all(.not.fail_max_list),&
  "Error find fields minmax equidistant mesh: maximum list not a maximum!")
  call assert_true(all(.not.fail_global),&
  "Error find fields minmax equisitant mesh: not global minimum or maximum!")
end subroutine test_field_minmax_equidistant_mesh

!> Tools ----------------------------------------------------
!>-----------------------------------------------------------
end module mod_fields_minmax_mpi_test
