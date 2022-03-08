!> mod_fields_minmax_test implements variables and procedures
!> used for testing methods use for finding jorek fields
!> minima and maxima
module mod_fields_minmax_mpi_test
use fruit
use constants,         only: TWOPI
use mod_rng,           only: type_rng
use data_structure,    only: type_node_list
use data_structure,    only: type_element_list
implicit none
private
public :: run_fruit_fields_minmax_mpi

!> Variables and datatypes ----------------------------------
character(len=31),parameter   :: filename_2d="test_jorek_2d_fields_restart.h5"
character(len=31),parameter   :: filename_3d="test_jorek_3d_fields_restart.h5"
integer,parameter             :: rst_format=0
integer,parameter             :: n_trials_sol=1000000
integer,parameter             :: n_trials_test=1000
integer,parameter             :: field_id_sol=1 !< poloidal flux
real*8,parameter              :: tol_real8=1.d-15
real*8,parameter              :: phi_sol=0.d0
real*8,dimension(2),parameter :: phi_interval_sol=(/0.d0,TWOPI/)
type(type_node_list)          :: node_list_2d,node_list_3d
type(type_element_list)       :: element_list_2d,element_list_3d
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
  call test_field_minmax_monte_carlo_2d(rank,n_tasks,ifail)
  call test_field_minmax_monte_carlo_3d(rank,n_tasks,ifail)
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
  call import_hdf5_restart(node_list_2d,element_list_2d,&
  trim(filename_2d),rst_format,ifail)
  call import_hdf5_restart(node_list_3d,element_list_3d,&
  trim(filename_3d),rst_format,ifail)
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
!> test the find of minimum and maximum values of a field 2d
subroutine test_field_minmax_monte_carlo_2d(rank,n_tasks,ifail)
  use mod_pcg32_rng,     only: pcg32_rng
  use mod_interp,        only: interp_PRZ
  use mod_fields_minmax, only: field_minmax_monte_carlo
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> variables:
  integer             :: ii,jj,kk
  real*8              :: R_test,Z_test
  real*8,dimension(1) :: field_test
  real*8,dimension(2) :: st,minmax_global
  real*8,dimension(2,element_list_2d%n_elements) :: minmax_list
  logical,dimension(element_list_2d%n_elements)  :: success_min_list
  logical,dimension(element_list_2d%n_elements)  :: success_max_list
  logical,dimension(2) :: success_global
  !> initialisation
  success_min_list = .true.; success_max_list = .true.; success_global = .true.;
  !> find critical points
  write(*,*) "searching for 2D field maxima and minima"
  call field_minmax_monte_carlo(node_list_2d,element_list_2d,&
  field_id_sol,n_trials_sol,(/phi_sol,phi_sol/),pcg32_rng(),rngs,minmax_list,&
  minmax_global,rank,n_tasks,ifail)
  !> checks if the extrema are the largest and smallest values for each point in the element
  write(*,*) "Checking 2D field maxima and minima"
  do ii=1,element_list_2d%n_elements
    do jj=1,n_trials_test
      call random_number(st)
      call interp_PRZ(node_list_2d,element_list_2d,ii,(/field_id_sol/),1,&
      st(1),st(2),phi_sol,field_test,R_test,Z_test)
      if(field_test(1).lt.minmax_list(1,ii)) success_min_list(ii) = .false.
      if(field_test(1).gt.minmax_list(2,ii)) success_max_list(ii) = .false.
      if(field_test(1).lt.minmax_global(1))  success_global(1) = .false.
      if(field_test(1).gt.minmax_global(2))  success_global(2) = .false.
    enddo
  enddo
  call assert_true(all(success_min_list),&
  "Error find 2D field minmax: minimum list not a minimum!")
  call assert_true(all(success_max_list),&
  "Error find 2D field minmax: maximum list not a maximum!")
  call assert_true(all(success_global),&
  "Error find 2D field minmax: not global minimum or maximum!")
end subroutine test_field_minmax_monte_carlo_2d

!> test the find of minimum and maximum values of a field 3d
subroutine test_field_minmax_monte_carlo_3d(rank,n_tasks,ifail)
  use mod_pcg32_rng,     only: pcg32_rng
  use mod_interp,        only: interp_PRZ
  use mod_fields_minmax, only: field_minmax_monte_carlo
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> variables:
  integer             :: ii,jj,kk
  real*8              :: R_test,Z_test
  real*8,dimension(1) :: field_test
  real*8,dimension(2) :: minmax_global
  real*8,dimension(3) :: stphi
  real*8,dimension(2,element_list_3d%n_elements) :: minmax_list
  logical,dimension(element_list_3d%n_elements)  :: success_min_list
  logical,dimension(element_list_3d%n_elements)  :: success_max_list
  logical,dimension(2) :: success_global
  !> initialisation
  success_min_list = .true.; success_max_list = .true.; success_global = .true.;
  !> find critical points
  write(*,*) "searching for 3D field maxima and minima"
  call field_minmax_monte_carlo(node_list_3d,element_list_3d,&
  field_id_sol,n_trials_sol,phi_interval_sol,pcg32_rng(),rngs,minmax_list,&
  minmax_global,rank,n_tasks,ifail)
  !> checks if the extrema are the largest and smallest values for each point in the element
  write(*,*) "Checking 3D field maxima and minima"
  do ii=1,element_list_3d%n_elements
    do jj=1,n_trials_test
      call random_number(stphi)
      stphi(3) = phi_interval_sol(1)  + (phi_interval_sol(2)-phi_interval_sol(1))*stphi(3)
      call interp_PRZ(node_list_3d,element_list_3d,ii,(/field_id_sol/),1,&
      stphi(1),stphi(2),stphi(3),field_test,R_test,Z_test)
      if(field_test(1).lt.minmax_list(1,ii)) success_min_list(ii) = .false.
      if(field_test(1).gt.minmax_list(2,ii)) success_max_list(ii) = .false.
      if(field_test(1).lt.minmax_global(1))  success_global(1) = .false.
      if(field_test(1).gt.minmax_global(2))  success_global(2) = .false.
    enddo
  enddo
  call assert_true(all(success_min_list),&
  "Error find 3D field minmax: minimum list not a minimum!")
  call assert_true(all(success_max_list),&
  "Error find 3D field minmax: maximum list not a maximum!")
  call assert_true(all(success_global),&
  "Error find 3D field minmax: not global minimum or maximum!")
end subroutine test_field_minmax_monte_carlo_3d

!> Tools ----------------------------------------------------
!>-----------------------------------------------------------
end module mod_fields_minmax_mpi_test
