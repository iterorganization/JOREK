!> particles/examples/proj_f.f90
!> Requires the JOREK input file
!>
!> CLI Arguments:
!>   * particle restart file name
!>   * (optional) jorek restart file name
!> Calculate projection of particles in elements of restart file name with
!> radiation
program proj_f_rad
use particle_tracer
use mod_project_particles
use gauss
use constants
use mod_interp
use mod_radiation, only: proj_Lz
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

allocate(particle_fieldline::sim%groups(1)%particles(n))
!!$omp parallel do default(none) private(i_p, i, j, k, l, R, R_s, R_t, Z, Z_s, Z_t, xjac) &
!!$    shared(sim)
do i_p=1,n_plane
  do i=1,sim%fields%element_list%n_elements
    do j=1,n_gauss
      do k=1,n_gauss
        l = (i_p-1) * n_gauss_2*sim%fields%element_list%n_elements + (i-1)*n_gauss_2 + (j-1)*n_gauss + k
        associate(p => sim%groups(1)%particles(l))
          p%x(3) = TWOPI*real(i_p,8)/real(n_plane,8)/real(n_period,8)
          p%i_elm = i
          p%st = [Xgauss(j), Xgauss(k)]
          ! Every particle represents a sample in the integral in an element
          ! with gaussian quadrature. We need the weights and the area here
          ! and a correction for the number of planes
          call interp_RZ(sim%fields%node_list,sim%fields%element_list,i,Xgauss(j), Xgauss(k), &
            R,R_s,R_t,Z,Z_s,Z_t)
          xjac = R_s*Z_t - R_t*Z_s
          p%weight = real(Wgauss(j)*Wgauss(k),4)*xjac*R*TWOPI/real(n_plane,4)
        end associate
      end do
    end do
  end do
end do
!!$omp end parallel do

! Set up the diagnostics output
proj = new_projection(sim%fields%node_list, sim%fields%element_list, filter=6d-5, &
    f=[proj_f(proj_Lz, group=1)], basename='qperp', &
    to_h5=.true.)
call with(sim, proj)

call sim%finalize
end program proj_f_rad
