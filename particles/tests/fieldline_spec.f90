!> This module contains some testcases for fieldline tracers
module fieldline_spec
use data_structure
use mod_particle_types
use fruit
use mod_sobseq_rng
use mod_initialise_particles
use mod_fields_linear
use projection_helpers
use mod_fieldline_euler
use mod_neighbours
implicit none

contains

!> Actions to perform before any of these tests
subroutine setup_fieldline_spec
  use basis_at_gaussian
  call initialise_basis !< Calculate the basis functions at the gaussian points
end subroutine setup_fieldline_spec

!> Test tracing a fieldline back and forth with euler
subroutine test_fieldline_backforth_euler
  use mod_find_rz_nearby
  type(jorek_fields_interp_linear) :: f
  integer, parameter :: n_p = 2
  type(particle_fieldline) :: p(n_p)
  real*8, parameter :: v = 1d5 ! 20 meters around the torus at this velocity
  real*8 :: rz_old(2), st_old(2), E(3), B(3), psi, U, psi0, dt, phi0
  real*8 :: R_out, Z_out, s_out, t_out ! against find_RZ trouble
  integer :: ielm_out ! against find_RZ trouble
  integer :: i, j, k, ifail, i_elm_old
  character(len=2) :: is
  if (.not. associated(f%node_list))    allocate(f%node_list)
  if (.not. associated(f%element_list)) allocate(f%element_list)
  call default_flux_grid_31(f%node_list, f%element_list)
  ! Call this once to setup the rtree
  call find_RZ(f%node_list,f%element_list,2.d0,1.d0,R_out,Z_out,ielm_out,s_out,t_out,ifail)
  ! Setup neighbour information for the run
  call update_neighbours(f%node_list, f%element_list)
  
  do i=-9,-7
    write(is,"(i2)") i
    dt = 10.d0**i
    call initialise_particles(p, f%node_list, f%element_list, sobseq_rng())
    do k=1,n_p
      p(k)%v = v
      call f%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi0, u)
      phi0 = p(k)%x(3)
      do j=1,20*10**(-(i+5)) ! 200 steps at smallest dt
        if (p(k)%i_elm .eq. 0) then
          call assert_true(.false., 'Particle should be in domain dt=10**'//is)
          exit
        end if
        call f%calc_EBpsiU(0.d0, p(k)%i_elm, &
            p(k)%st, p(k)%x(3), E, B, psi, U)
        rz_old    = p(k)%x(1:2)
        st_old    = p(k)%st
        i_elm_old = p(k)%i_elm

        call fieldline_euler_push_cylindrical(p(k), B, dt)
        call find_RZ_nearby(f%node_list, f%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
            p(k)%x(1), p(k)%x(2), p(k)%st(1), p(k)%st(2), p(k)%i_elm, ifail)
        if (j .eq. 10*10**(-(i+5))) then
          call assert_equals(0.d0, psi-psi0, 13d4*dt, "Must not leave flux surface mid dt=1e"//is)
          p(k)%v = -v ! go backwards after this point
        end if
      end do
      if (p(k)%i_elm .ne. 0) then
        call f%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi, u)
        call assert_equals(0.d0, psi-psi0, 25d4*dt, "Must not leave flux surface dt=1e"//is)
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
  type(jorek_fields_interp_linear) :: f
  integer, parameter :: n_p = 2
  type(particle_fieldline) :: p(n_p)
  real*8, parameter :: v = 1d5 ! 20 meters around the torus at this velocity
  real*8 :: rz_old(2), st_old(2), E(3), B(3), psi, U, psi0, dt, phi0
  real*8 :: R_out, Z_out, s_out, t_out ! against find_RZ trouble
  integer :: ielm_out ! against find_RZ trouble
  integer :: i, j, k, ifail, i_elm_old
  character(len=2) :: is
  if (.not. associated(f%node_list))    allocate(f%node_list)
  if (.not. associated(f%element_list)) allocate(f%element_list)
  call default_flux_grid_31(f%node_list, f%element_list)
  ! Call this once to setup the rtree
  call find_RZ(f%node_list,f%element_list,2.d0,1.d0,R_out,Z_out,ielm_out,s_out,t_out,ifail)
  ! Setup neighbour information for the run
  call update_neighbours(f%node_list, f%element_list)
  
  do i=-9,-6
    write(is,"(i2)") i
    dt = 10.d0**i
    call initialise_particles(p, f%node_list, f%element_list, sobseq_rng())
    do k=1,n_p
      p(k)%v = v
      call f%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi0, u)
      phi0 = p(k)%x(3)
      ! Do a single euler step forward to setup the adams-bashforth method
      rz_old    = p(k)%x(1:2)
      st_old    = p(k)%st
      i_elm_old = p(k)%i_elm
      call fieldline_euler_push_cylindrical(p(k), B, dt)
      call find_RZ_nearby(f%node_list, f%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
          p(k)%x(1), p(k)%x(2), p(k)%st(1), p(k)%st(2), p(k)%i_elm, ifail)
      p(k)%B_hat_prev = B/norm2(B)

      do j=2,20*10**(-(i+5)) ! 200 steps at smallest dt
        if (p(k)%i_elm .eq. 0) then
          call assert_true(.false., 'Particle should be in domain dt=1e'//is)
          exit
        end if
        call f%calc_EBpsiU(0.d0, p(k)%i_elm, &
            p(k)%st, p(k)%x(3), E, B, psi, U)
        rz_old    = p(k)%x(1:2)
        st_old    = p(k)%st
        i_elm_old = p(k)%i_elm

        call fieldline_adams_bashforth_push_cylindrical(p(k), B, dt)
        call find_RZ_nearby(f%node_list, f%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
            p(k)%x(1), p(k)%x(2), p(k)%st(1), p(k)%st(2), p(k)%i_elm, ifail)
        if (j .eq. 10*10**(-(i+5))) then
          call assert_equals(0.d0, psi-psi0, 2d9*dt**2, "Must not leave flux surface mid dt=1e"//is)
          p(k)%v = -v ! go backwards after this point
        end if
      end do
      if (p(k)%i_elm .ne. 0) then
        call f%calc_EBpsiU(0.d0, p(k)%i_elm, p(k)%st, p(k)%x(3), E, B, psi, u)
        call assert_equals(0.d0, psi-psi0, 2d9*dt**2, "Must not leave flux surface dt=1e"//is)
        call assert_equals(0.d0, p(k)%x(3)-phi0, 8d4*dt, "Must be back at same phi dt=1e"//is) ! WARNING: this is linear instead of quadratic
      else
        call assert_true(.false., 'Particle should be in domain after run dt=1e'//is)
      end if
    end do
  end do
end subroutine test_fieldline_backforth_adams_bashforth

end module fieldline_spec
