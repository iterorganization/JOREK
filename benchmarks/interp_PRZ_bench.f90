!> Benchmark program to test interpolation speeds of mod_interp_PRZ
!> Set n to switch element and position every n steps
!>
!> run with ulimit -s unlimited
program interp_PRZ_bench
use data_structure
use mpi
use mod_pcg32_rng
use mod_random_seed
use mod_interp
use mod_parameters, only: n_degrees
use basis_at_gaussian, only: initialise_basis
implicit none

integer, parameter :: n_switch = 100 !< Go to a new element every 100 steps (mimic jorek particle pattern)
! we go to a new position inside the element every step.
integer, parameter :: N_loops = 1000000
integer, parameter :: n_v = 1, i_var(n_v) = [1]

integer :: ierr, provided, i, j, k, i_elm
type(type_node_list) :: node_list
type(type_element_list) :: element_list
real*8 :: t0, t1
type(pcg32_rng) :: rng
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
real*8, dimension(n_v) :: P, P_s, P_t, P_phi
real*8 :: u(1), s, t, phi

real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
integer, parameter :: n = 100 ! number of nodes in r, z directions

! Setup the grid
call MPI_INIT_THREAD(MPI_THREAD_SINGLE, provided, ierr)
call initialise_basis
node_list%n_nodes = 0
element_list%n_elements = 0
call grid_bezier_square(n, n, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)

! Set random values here
call rng%initialize(1, random_seed(), 1, 1)

do i=1,node_list%n_nodes
  do k=1,n_degrees
    call rng%next(u)
    node_list%node(i)%values(1,k,:) = u(1)
    node_list%node(i)%deltas(1,k,:) = 0.1d0*u(1)
  end do
enddo

! Wait after this omp part for a while
call sleep(3)

call cpu_time(t0)
! Start the test
i_elm = 1
do i=1,N_loops
  if (mod(i,n_switch) == 0 .and. n_switch > 0) then
    call rng%next(u)
    i_elm = nint(1.d0+u(1)*real(element_list%n_elements-1))
  end if
  call rng%next(u); s = u(1)
  call rng%next(u); t = u(1)
  call rng%next(u); phi = u(1)
  !call interp_PRZ(node_list, element_list, i_elm, i_var, n_v, s, t, phi, P, P_s, &
  !P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t)
  call interp_RZ(node_list, element_list, i_elm, s, t, R, R_s, R_t, Z, Z_s, Z_t)
  call rng%next(u)
  do k=1,n_v
    if (u(1) .gt. 0.999999) write(*,"(10g16.7)") P(k), P_s(k), P_t(k), P_phi(k), R, R_s, R_t, Z, Z_s, Z_t
  end do
end do
call cpu_time(t1)

write(*,*) "Time for ", N_loops, " iterations: ", t1-t0, "s"
write(*,*) "Time per iteration: ", ((t1-t0)/real(N_loops)) *1d6, " microseconds"

end program interp_PRZ_bench
