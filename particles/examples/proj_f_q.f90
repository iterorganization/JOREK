!> particles/examples/proj_f.f90
!> Requires: jorek_restart.h5 and the JOREK input file
!> Read fields from jorek_restart.h5 and calculate projection of some function
!> evaluated at gauss points.
program proj_f_q
use particle_tracer
use mod_project_particles
use gauss
use constants
use mod_interp
!$ use omp_lib
implicit none

type(projection) :: proj

type(event) :: fieldreader
integer :: i, j, k, l, i_p, tid, n_threads, n
character(len=20) :: time_s
real*8 :: R, R_s, R_t, Z, Z_s, Z_t, xjac

! Start up MPI, jorek
call sim%initialize(num_groups=1)

!call get_command_argument(1, time_s)
!read(time_s,*) sim%time

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(i=-1))
call with(sim, fieldreader)

! Set up particles
sim%groups(:)%Z    = -1
sim%groups(:)%mass = -1

n = sim%fields%element_list%n_elements * n_gauss_2 * n_plane

!!$omp parallel do default(none) private(i_p, i, j, k, l, R, R_s, R_t, Z, Z_s, Z_t, xjac) &
!!$    shared(sim)
do i_p=1,n_plane
  do i=1,sim%fields%element_list%n_elements
    do j=1,n_gauss
      do k=1,n_gauss
        l = (i_p-1) * n_gauss_2*sim%fields%element_list%n_elements + (i-1)*n_gauss_2 + (j-1)*n_gauss + k
        sim%groups(1)%particles(l)%x(3) = TWOPI*real(i_p,8)/real(n_plane,8)/real(n_period,8)
        sim%groups(1)%particles(l)%i_elm = i
        sim%groups(1)%particles(l)%st = [Xgauss(j), Xgauss(k)]
        ! Every particle represents a sample in the integral in an element
        ! with gaussian quadrature. We need the weights and the area here
        ! and a correction for the number of planes
        call interp_RZ(sim%fields%node_list,sim%fields%element_list,i,Xgauss(j), Xgauss(k), &
          R,R_s,R_t,Z,Z_s,Z_t)
        xjac = R_s*Z_t - R_t*Z_s
        sim%groups(1)%particles(l)%weight = real(Wgauss(j)*Wgauss(k),4)*xjac*R*TWOPI/real(n_plane,4)
      end do
    end do
  end do
end do
!!$omp end parallel do

! Set up the diagnostics output
proj = new_projection(sim%fields%node_list, sim%fields%element_list, smoothing=6d-5, &
    f=[proj_f(proj_q_homma, 1)], basename='qperp', &
    to_h5=.true.)
call with(sim, proj)

call sim%finalize
contains
!> Interpolate the parallel electric field at the particle position
pure function proj_q_homma(sim, group, particle)
  use constants
  use mod_collisions
  use mod_ionisation_recombination
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_q_homma, E(3), B(3), psi, U, q(3)
  real*8 :: n_e, T_e, grad_T_e(3) ! temperature and gradient in [K]
  call sim%fields%calc_NeTe(sim%time, particle%i_elm, particle%st, particle%x(3),n_e, T_e, grad_T_e)
  call sim%fields%calc_EBpsiU(sim%time, particle%i_elm, particle%st, particle%x(3), &
      E, B, psi, U)
  q = q_homma2013(T_e*K_BOLTZ, grad_T_e*K_BOLTZ, B, n_e, 2.5d0, 1_1)

  ! choose your component
  ! hardcode first component
  proj_q_homma = q(1) !norm2(q - B*dot_product(q,B)/dot_product(B,B)) ! perp, total
  !proj_q = dot_product(q,B)/norm2(B) ! par
end function proj_q_homma
end program Proj_f_q
