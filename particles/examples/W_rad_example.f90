!> particles/examples/W_rad_example.f90
!> Requires: acd50_w.dat, scd50_w.dat plt50_w.dat, prb50_w.dat, jorek_restart.h5 and the JOREK input file
!> Compile with `j2p W_rad_example` and run with `W_rad_example < in1`
!> Read fields from jorek_restart.h5 and calculate radiated power
program W_rad_example
use particle_tracer
use mod_sobseq_rng
use mod_particle_diagnostics
use mod_project_particles
use mod_radiation
use mod_particle_types
implicit none

type(coronal)   :: cor
type(projection) :: proj

real*8, parameter :: alpha = 0.d0 !< if 0, sample from maxwellian
type(event) :: fieldreader
integer :: i
integer :: ierr, my_id
real*8, parameter :: timesteps(1) = [1d-9]

! Start up MPI, jorek
call sim%initialize(num_groups=1)

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(&
    basename='jorek', i=-1))
call with(sim, fieldreader)

! Prepare the coronal equilibrium
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
sim%groups(1)%ad = read_adf11(my_id, '50_w')
cor = coronal(sim%groups(1)%ad)

! Set up particles
sim%groups(:)%Z    = 74
sim%groups(:)%mass = 183.84 !< atomic mass units
do i=1,1
  allocate(particle_kinetic_leapfrog::sim%groups(i)%particles(10000000))
  ! For every particle accept or reject it with probability f_psi_inside(psi)
  call initialise_particles_H_mu_psi(sim%groups(i)%particles, &
      sim%fields, sobseq_rng(), &
      sim%groups(i)%mass, &
      uniform_space=.true., cor=cor)

  select type (p => sim%groups(i)%particles)
  type is (particle_kinetic_leapfrog)
    p(:)%weight = 1.d0
    call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(i)%mass, &
        sim%fields, sim%time, timesteps(i))
  end select
end do

call adjust_particle_weights(sim%groups(1)%particles, 6d16) ! 60 m^3 at n_i = 10^15 /m^3

! Set up the diagnostics output
proj = projection(sim%fields%node_list, sim%fields%element_list, filter=1d-5, &
    f=[proj_f(proj_Lz, 1)], &
    to_h5=.true., to_vtk=.true.)

call with(sim, proj)

call sim%finalize
end program W_rad_example
