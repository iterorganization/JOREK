!> Given a jorek grid initialize particles and write out a particle restart file.
!> set a fixed charge.
!> needs jorek_restart.h5 and run with < input file
program seed_particles
use particle_tracer
use mod_particle_io
use mod_rej_f
use initialisers_base
implicit none

type(event) :: fieldreader
type(spatial_pdf) :: space_pdf
integer :: i, j

! Start up MPI, jorek
call sim%initialize(num_groups=1)

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)

! Set up spatial pdf type
space_pdf%f    => f_psi_inside
space_pdf%vars =  [1]

! Set up particles
sim%groups(:)%Z    = -2 
sim%groups(:)%mass = 2.d0 !< atomic mass units
do i=1,1
  allocate(particle_kinetic_leapfrog::sim%groups(i)%particles(2000))
  ! For every particle accept or reject it with probability f_psi_inside(psi)
  ! this function below overwrites the whole particle set
  call initialise_particles_H_mu_psi(sim%groups(i)%particles, &
      sim%fields, pcg32_rng(), &
      sim%groups(i)%mass, &
      uniform_space=.true., space_pdf=space_pdf, charge=1)
  sim%groups(i)%particles(:)%weight = 1.0
end do

call write_simulation_hdf5(sim, 'part_restart.h5')
contains
pure function f_psi_inside(n, P, grad_P) result(f)
  integer,                intent(in) :: n
  real*8, dimension(n),   intent(in) :: P
  real*8, dimension(3,n), intent(in) :: grad_P
  real*4 :: f
  f = 5e-3
  if (P(1) .lt. -0.26 .and. P(1) .gt. -0.35) f = 1e0
end function f_psi_inside
end program seed_particles
