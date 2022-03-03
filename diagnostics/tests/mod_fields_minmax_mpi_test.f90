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
  call test_generation_equidistant_stphi_mesh(rank,n_tasks,ifail)
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
subroutine test_generation_equidistant_stphi_mesh(rank,n_tasks,ifail)
  use constants,         only: TWOPI
  use mod_fields_minmax, only: generate_stphi_mesh
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  !> variables
  integer :: ii
  real*8 :: d_s_sol, d_t_sol,d_phi_sol
  real*8,dimension(3) :: d_test
  real*8,dimension(3,n_s*n_t*n_phi)    :: eq_mesh
  logical,dimension(3,n_s*n_t*n_phi-1) :: success
  !> initialisations
  d_s_sol = 1.d0/real(n_s-1,kind=8); d_t_sol = 1.d0/real(n_t-1,kind=8);
  d_phi_sol = TWOPI/real(n_phi-1,kind=8); success = .true.;
  !> generate equidistant mesh
  call generate_stphi_mesh(n_s,n_t,n_phi,eq_mesh)
  !> checks
  do ii=1,n_s*n_t*n_phi-1
    d_test = eq_mesh(:,ii+1)-eq_mesh(:,ii)
    if(d_test(1).eq.-1.d0) d_test(1) = 0.d0
    if(d_test(2).eq.-1.d0) d_test(2) = 0.d0
    if(d_test(1).ne.0.d0) success(1,ii) = (d_test(1).ge.(d_s_sol-tol_real8)).and.&
                                          (d_test(1).le.(d_s_sol+tol_real8))
    if(d_test(2).ne.0.d0) success(2,ii) = (d_test(2).ge.(d_t_sol-tol_real8)).and.&
                                          (d_test(2).le.(d_t_sol+tol_real8))
    if(d_test(3).ne.0.d0) success(3,ii) = (d_test(3).ge.(d_phi_sol-tol_real8)).and.&
                                          (d_test(3).le.(d_phi_sol+tol_real8))
  enddo
  call assert_true(all(success),"Error fields minmax equidistant mesh: mesh size mismatch!")
end subroutine test_generation_equidistant_stphi_mesh

!> test the generation of random mesh
!> test the find of extrema on equidistant mesh
!> test the find of extrema on random mesh
!> Tools ----------------------------------------------------
!>-----------------------------------------------------------
end module mod_fields_minmax_mpi_test
