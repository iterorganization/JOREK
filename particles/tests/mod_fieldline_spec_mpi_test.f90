!> This module contains some testcases for fieldline tracers
module mod_fieldline_spec_mpi_test
use fruit
use fruit_mpi
use data_structure
use mod_particle_types
use mod_sobseq_rng
use mod_initialise_particles
use initialisers_base
use mod_fields_linear
use mod_fieldline_euler
use mod_neighbours
implicit none
private
public :: run_fruit_fieldline_spec_mpi
!> Variables --------------------------------------
integer,parameter :: master_rank=0
integer,parameter :: n_poloidal_nodes=31
integer,parameter :: n_radial_nodes=30
real*8,parameter  :: R_init=1.5d0
real*8,parameter  :: Z_init=2.d-1
integer :: rank_loc,n_tasks_loc,ifail_loc
type(type_node_list),target      :: node_list_sol
type(type_element_list),target   :: element_list_sol
type(jorek_fields_interp_linear) :: fields_sol
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_fieldline_spec_mpi(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in) :: rank,n_tasks
  if(rank.eq.master_rank) write(*,'(/A)') "  ... setting-up: fieldline spec"
  if(rank.eq.master_rank) write(*,'(/A)') "  ... running: fieldline spec"
  call setup(rank,n_tasks,ifail)
  call run_test_case(test_fieldline_backforth_euler,'test_fieldline_backforth_euler')
  call run_test_case(test_fieldline_backforth_adams_bashforth,'test_fieldline_backforth_adams_bashforth')
  if(rank.eq.master_rank) write(*,'(/A)') "  ... tearing-down: fieldline spec"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_fieldline_spec_mpi

!> Set-up and tear-down ---------------------------
!> Actions to perform before any of these tests
subroutine setup(rank,n_tasks,ifail)
  use mod_projection_helpers_test_tools, only: default_flux_grid
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc = rank; n_tasks_loc = n_tasks; ifail_loc = ifail;
  call init_node_list(node_list_sol, n_nodes_max, node_list_sol%n_dof, n_var)

  fields_sol%node_list    => node_list_sol
  fields_sol%element_list => element_list_sol
  !> compute the MHD equilibrium and define the flux grid
  call default_flux_grid(rank_loc,n_tasks_loc,n_poloidal_nodes,&
  n_radial_nodes,fields_sol%node_list,fields_sol%element_list,ifail_loc)
end subroutine setup

subroutine teardown(rank,n_tasks,ifail)
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc = 0; n_tasks_loc = 0; ifail = ifail_loc;
  if(associated(fields_sol%node_list))    fields_sol%node_list => NULL()
  if(associated(fields_sol%element_list)) fields_sol%element_list => NULL()
  call dealloc_node_list(node_list_sol)
end subroutine teardown

!> Tests ------------------------------------------
!> Test tracing a fieldline back and forth with euler
subroutine test_fieldline_backforth_euler
  use mod_find_rz_nearby
  integer, parameter       :: n_p = 2
  type(particle_fieldline) :: p(n_p)
  real*8, parameter        :: v = 1d5 ! 20 meters around the torus at this velocity
  real*8                   :: rz_old(2), st_old(2), E(3), B(3), psi, U, psi0, dt, phi0
  real*8                   :: R_out, Z_out, s_out, t_out ! against find_RZ trouble
  integer                  :: ielm_out ! against find_RZ trouble
  integer                  :: i, j, k, ifail, i_elm_old
  character(len=2)         :: is
  
  ! Call this once to setup the rtree
  call find_RZ(fields_sol%node_list,fields_sol%element_list,R_init,Z_init,&
  R_out,Z_out,ielm_out,s_out,t_out,ifail)
  
  ! Setup neighbour information for the run
  call update_neighbours(fields_sol%node_list, fields_sol%element_list)
  
  do i=-9,-7
    
    write(is,"(i2)") i
    dt = 10.d0**i
    call initialise_particles(p, fields_sol%node_list, fields_sol%element_list, sobseq_rng())
    do k=1,n_p
    
      p(k)%v = v
      
      call fields_sol%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi0, u)
      
      phi0 = p(k)%x(3)
      
      do j=1,20*10**(-(i+5)) ! 200 steps at smallest dt

        if (p(k)%i_elm .le. 0) then
          call assert_true(.false., 'Particle should be in domain dt=10**'//is)
          exit
        end if

        rz_old    = p(k)%x(1:2)
        st_old    = p(k)%st
        i_elm_old = p(k)%i_elm

        call fieldline_euler_push_cylindrical(p(k), B, dt)

        call find_RZ_nearby(fields_sol%node_list, fields_sol%element_list, &
             rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
            p(k)%x(1), p(k)%x(2), p(k)%st(1), p(k)%st(2), p(k)%i_elm, ifail)

        call fields_sol%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi, U)

        if (j .eq. 10*10**(-(i+5))) then
          call assert_equals(0.d0, psi-psi0, 16d4*dt, "Must not leave flux surface mid dt=1e"//is)
          p(k)%v = -v ! go backwards after this point
        end if

      end do

      if (p(k)%i_elm .gt. 0) then
        call fields_sol%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi, u)
        call assert_equals(0.d0, psi-psi0, 32d4*dt, "Must not leave flux surface dt=1e"//is)
        call assert_equals(0.d0, p(k)%x(3)-phi0, 7d6*dt, "Must be back at same phi dt=1e"//is)
      else
        call assert_true(.false., 'Particle should be in domain after run dt=1e'//is)
      end if

    end do
  end do
end subroutine test_fieldline_backforth_euler

!> Test tracing a fieldline back and forth with Adams-Bashforth
subroutine test_fieldline_backforth_adams_bashforth
  use mod_find_rz_nearby
  integer, parameter       :: n_p = 2
  type(particle_fieldline) :: p(n_p)
  real*8, parameter        :: v = 1d5 ! 20 meters around the torus at this velocity
  real*8                   :: rz_old(2), st_old(2), E(3), B(3), psi, U, psi0, dt, phi0
  real*8                   :: R_out, Z_out, s_out, t_out ! against find_RZ trouble
  integer                  :: ielm_out ! against find_RZ trouble
  integer                  :: i, j, k, ifail, i_elm_old
  character(len=2)         :: is
 
  ! Call this once to setup the rtree
  call find_RZ(fields_sol%node_list,fields_sol%element_list,R_init,Z_init,&
   R_out,Z_out,ielm_out,s_out,t_out,ifail)
  
  ! Setup neighbour information for the run
  call update_neighbours(fields_sol%node_list, fields_sol%element_list)
  
  do i=-9,-6
    
    write(is,"(i2)") i
    dt = 10.d0**i
    
    call initialise_particles(p, fields_sol%node_list, fields_sol%element_list, sobseq_rng())

    do k=1,n_p
    
      p(k)%v = v
      
      call fields_sol%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi0, u)
      
      phi0 = p(k)%x(3)
      ! Do a single euler step forward to setup the adams-bashforth method
      rz_old    = p(k)%x(1:2)
      st_old    = p(k)%st
      i_elm_old = p(k)%i_elm

      call fieldline_euler_push_cylindrical(p(k), B, dt)
      
      call find_RZ_nearby(fields_sol%node_list, fields_sol%element_list, &
           rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
           p(k)%x(1), p(k)%x(2), p(k)%st(1), p(k)%st(2), p(k)%i_elm, ifail)
      
      p(k)%B_hat_prev = B/norm2(B)

      do j=2,20*10**(-(i+5)) ! 200 steps at smallest dt
      
        if (p(k)%i_elm .le. 0) then
          call assert_true(.false., 'Particle should be in domain dt=1e'//is)
          exit
        end if

        call fields_sol%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi, U)

        rz_old    = p(k)%x(1:2)
        st_old    = p(k)%st
        i_elm_old = p(k)%i_elm

        call fieldline_adams_bashforth_push_cylindrical(p(k), B, dt)
        
        call find_RZ_nearby(fields_sol%node_list, fields_sol%element_list, &
             rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
             p(k)%x(1), p(k)%x(2), p(k)%st(1), p(k)%st(2), p(k)%i_elm, ifail)
        
        if (j .eq. 10*10**(-(i+5))) then
          call assert_equals(0.d0, psi-psi0, 2d9*dt**2, "Must not leave flux surface mid dt=1e"//is)
          p(k)%v = -v ! go backwards after this point
        end if
      
      end do
      
      if (p(k)%i_elm .gt. 0) then
        call fields_sol%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi, u)
        call assert_equals(0.d0, psi-psi0, 2d9*dt**2, "Must not leave flux surface dt=1e"//is)
        call assert_equals(0.d0, p(k)%x(3)-phi0, 8d4*dt, "Must be back at same phi dt=1e"//is) ! WARNING: this is linear instead of quadratic
      else
        call assert_true(.false., 'Particle should be in domain after run dt=1e'//is)
      end if

    end do
  end do
end subroutine test_fieldline_backforth_adams_bashforth

!> ------------------------------------------------
end module mod_fieldline_spec_mpi_test
