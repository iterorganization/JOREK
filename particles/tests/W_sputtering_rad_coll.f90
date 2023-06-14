!> Test program to sputter particles from a number of point sources
!> with static fields, with collisions, projection and radiation diagnostics
!>
!> We make many simplifying assumptions here
!> - collision dt is fixed to 10^-10, should depend on plasma params?
!> - particle sources are a series of points
!> - We do 10 collisions from each maxwellian to save JOREK field evaluation costs. Assume this is ok
program W_sputtering_rad_coll
use particle_tracer
use mod_sobseq_rng
use mod_pcg32_rng
use mod_random_seed
use mod_particle_diagnostics
use mod_project_particles
use mod_ionisation_recombination
use mpi
use mod_sampling
use mod_particle_io
use mod_radiation, only: proj_Lz
use mod_sputtering, only: simple_sputter !< simple sputtering implementation
use phys_module, only: central_density
use mod_collisions
use mod_interp, only: interp_RZ
!$ use omp_lib
implicit none

type(adf11_all) :: adas
integer, parameter :: n_sources = 20
! Take n_sources positions, from -4.22949 down to -4.51725
real*8, parameter :: top = -4.22949d0, bottom = -4.51725d0
real*8, parameter, dimension(n_sources) :: YiCi = [&
    0.d0, 0.d0, 0.d0, 0.d0, 1d-5, 7d-5, 2d-4, 4d-4, 8d-4, 1d-3, &
    1.5d-3, 1.9d-3, 2.3d-3, 3d-3, 3.2d-3, 2.8d-3, 9d-4, 1d-5, 0.d0, 0.d0]
type(simple_sputter), dimension(n_sources) :: sources

integer, parameter :: n_coll = 100 ! number of collisions at any specific point to calculate (reference point: 1 per 10^-10 s)
integer(kind=1), parameter :: q_b = 1_1
real*8, parameter :: m_b = 2.d0 !< [u]

real*8 :: x(2), s, R, Z, R_s, R_t, Z_s, Z_t, c
real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
integer :: i_elm

real*8 :: timesteps(1) = [1d-8]
integer :: i, j, k, l, n_steps, n_lost, n_lost_all, i_elm_old, ifail, n_stream, seed, i_rng, ierr
real*8 :: target_time, t, n_e, T_e, grad_T_e(3), ran(6), q(3), coulomb_log, kTb, n_b, v_b(3,n_coll), ran2(6,n_coll)
!$ real*8 :: w0, w1
real*8 :: E(3), B(3), rz_old(2), st_old(2), psi, U, mmm(3)

type(event) :: fieldreader
type(pcg32_rng), dimension(:), allocatable :: rng

logical, parameter :: restart = .false.

! Start up MPI, jorek
call sim%initialize(num_groups=1)

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)

! Prepare the coronal equilibrium
adas = read_adf11(sim%my_id,'50_w')

! Set up particles
if (.not. restart) then
  allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(1000000))
  do i=1,size(sim%groups(1)%particles,1)
    sim%groups(1)%particles(i)%i_elm = 0
  end do
else
  call read_simulation_hdf5(sim, 'part_restart.h5')
end if
sim%groups(1)%Z    = 74
sim%groups(1)%mass = 183.84 !< atomic mass units
  

! Set up sputtering source at center of element containing position below
do i=1,n_sources
  x = [5.56121d0,top + (bottom - top)*real(i,8)/real(n_sources,8)]
  ! Find the closest position there
  call find_RZ(sim%fields%node_list, sim%fields%element_list, x(1), x(2), R, Z, &
      i_elm, s, t, ifail)
  ! Manually lookup sputtering coefficients with temperature given below (I know, I'm sorry...)
  ! It would be much much better to have some guess of the plasma composition and use that to derive the rates
  call sim%fields%calc_NeTe(sim%time, i_elm, [s,t], 0.d0, n_e, T_e) ! n_e in m^-3, T_e in K
  ! See Eckstein, 2008, Vacuum
  ! Note that we measure these values slightly away from the edge, as we need the temperature at the pre-sheath.

  ! Set t to 0 to go exactly to the edge for the sputtering coordinates
  t = 0.d0
  ! Recalculate R, Z
  call interp_RZ(sim%fields%node_list, sim%fields%element_list, i_elm, s, t, R, R_s, R_t, Z, Z_s, Z_t)

  !write(*,*) i, T_e*K_BOLTZ/EL_CHG ! [eV]

  ! The sputtering source is given by c*n*dx
  ! Set sputtering rate
  ! width ~ 25 cm, density ~ 0.12 * n0
  ! T ~ 0.003 * (Kb * mu_0 * n0)
  sources(i) = simple_sputter(x=[R,Z], st=[s,t], i_elm=i_elm, normal=[-1.d0, 0.d0, 0.d0], &
      source=0.25d0/real(n_sources)*YiCi(i)*n_e, &
      sim_source=3d7*YiCi(i)/maxval(YiCi), &
      ! calculate this from the YiCi to make particles with roughly equal weight
      ! make 0.01 billion particles per second in this source
      ! (*1e-6 = 20 per sputtering event per source per cpu at most)
      fill_fraction=0.1d0, &
      E_dist=thompson_dist(E_b=8.7d0, n=2), &
      last_time=sim%time)
end do

    

! Add only sources where there are nonzero c, i.e. 3..8 inclusive
events = [fieldreader, &
          event(sources(5), step=1d-6), &
          event(sources(6), step=1d-6), &
          event(sources(7), step=1d-6), &
          event(sources(8), step=1d-6), &
          event(sources(9), step=1d-6), &
          event(sources(10), step=1d-6), &
          event(sources(11), step=1d-6), &
          event(sources(12), step=1d-6), &
          event(sources(13), step=1d-6), &
          event(sources(14), step=1d-6), &
          event(sources(15), step=1d-6), &
          event(sources(16), step=1d-6), &
          event(sources(17), step=1d-6), &
          event(sources(18), step=1d-6), &
          event(write_particle_diagnostics(filename='diag.h5', append=restart), step=1d-5), &
          event(write_action(),   step=1d-5), &
          !event(projection(sim%fields%node_list, sim%fields%element_list, &
          !  filter=5d-6, proj_f=proj_Lz, to_h5=.true.),  step=1d-5), &
          event(projection(sim%fields%node_list, sim%fields%element_list, &
             filter=5d-6, to_h5=.true.),  step=1d-5) &
        ]
!call check_and_fix_timesteps(timesteps, events)

! Setup RNG for many threads
seed = random_seed()
n_stream = 1
!$ n_stream = omp_get_max_threads()
allocate(rng(n_stream))
do i=1,n_stream
  call rng(i)%initialize(6, seed, n_stream, i)
end do

do while (.not. sim%stop_now)
  target_time = next_event_at(sim, events)

  !$ w0 = omp_get_wtime()
  do i=1,size(sim%groups,1)
    n_steps = floor((target_time - sim%time)/timesteps(i))
    n_lost = 0

    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      !$omp parallel default(private), shared(sim, n_steps, timesteps, i, rng, adas, central_density)
      i_rng = 1
      !$ i_rng = omp_get_thread_num()+1

      !$omp do schedule(dynamic)
      do j=1,size(particles,1)
        do k=1,n_steps
          if (particles(j)%i_elm .le. 0) exit
          t = sim%time + k*timesteps(i)
          call sim%fields%calc_EBpsiU(t, particles(j)%i_elm, &
              particles(j)%st, particles(j)%x(3), E, B, psi, U)
          rz_old    = particles(j)%x(1:2)
          st_old    = particles(j)%st
          i_elm_old = particles(j)%i_elm

          ! Update charge based on ionisation coefficients
          call sim%fields%calc_NeTe(t, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), n_e, T_e, grad_T_e)
          call rng(i_rng)%next(ran)
          particles(j)%q = int(new_charge(int(particles(j)%q,4), adas, log10(n_e), log10(T_e), timesteps(i), ran(1:2)),1)

          if (particles(j)%q .gt. 0) then
            ! Calculate collisions
            n_b = n_e
            q = q_homma2013(T_e*K_BOLTZ, grad_T_e*K_BOLTZ, B, n_b, m_b, q_b)
            coulomb_log = coulomb_logarithm(kTb, n_b, particles(j)%q, q_b, sim%groups(i)%mass, m_b)
            ! Get parallel flow velocity
            call sim%fields%interp_PRZ(t, particles(j)%i_elm, [7], 1, particles(j)%st(1), particles(j)%st(2), &
                particles(j)%x(3), P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

            do l=1,n_coll
              call rng(i_rng)%next(ran2(:,l))
            end do
            call sample_velocity_dist_magnetized(n_coll, ran2(1:6,:), kTb, q, n_b, m_b, q_b, P(1)*B/sim%t_norm, v_b)

            do l=1,n_coll
              call rng(i_rng)%next(ran)
              call collide_particles(ran(1:3), particles(j)%q, sim%groups(i)%mass, particles(j)%v, &
                  q_b, m_b, v_b(:,l), n_b, coulomb_log, timesteps(i)/real(n_coll,8))
            end do
          end if

          call boris_push_cylindrical(particles(j), sim%groups(i)%mass, E, B, timesteps(i))
          call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
              particles(j)%x(1), particles(j)%x(2), particles(j)%st(1), particles(j)%st(2), particles(j)%i_elm, ifail)
        end do ! steps
      end do ! particles
      !$omp end do
      !$omp end parallel
    end select
  end do ! groups
  !$ w1 = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w1-w0)
  !$ if (sim%my_id .eq. 0) write(*,"(f10.7,A,3f9.4,A)") sim%time, " Particle stepping complete in ", mmm, "s"
  sim%time = target_time
  call with(sim, events, at=sim%time)
end do
call sim%finalize

end program W_sputtering_rad_coll
