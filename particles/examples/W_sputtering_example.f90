program W_sputtering_example

!> Test program to sputter particles from the divertor
!> with static or dynamic fields, with collisions, projection and radiation diagnostics
!>
!> We make many simplifying assumptions here
!> - collision dt is fixed to 10^-10, should depend on plasma params
!> - We do 10 collisions from each maxwellian to save JOREK field evaluation costs. Assume this is ok
use particle_tracer
use mod_sobseq_rng
use mod_pcg32_rng
use mod_random_seed
use mod_particle_diagnostics
use mod_project_particles
use mod_ionisation_recombination
use mpi
use mod_sampling
use mod_radiation, only: proj_Lz, proj_Lz_equil
use mod_particle_sputtering, only: particle_sputter
use phys_module, only: central_density
use mod_collisions
use mod_coronal
use mod_atomic_elements
use mod_particle_io

use mod_edge_domain
use mod_edge_elements
!$ use omp_lib
implicit none

type(particle_sputter) :: W_sputter_source_in, W_sputter_source_out

integer, parameter :: n_coll = 100 ! number of collisions at any specific point to calculate (reference point: approx 1 per 10^-10 s)
integer, parameter :: n_particles = 4000000
real*8, parameter :: sputter_fill_time = 2d-5 ! Time for the particle inventory to be full from the sputtering yield given no recycling
real*8, parameter :: sputter_dt = 1d-6 !< Time between sputtering actions
real*8 :: timesteps(2) = [1d-8,1d-8] !<1d-8
real*8 :: step_rest_time(2) = 0.d0 !< [s]
integer(kind=1), parameter :: q_b = 1_1
real*8, parameter :: m_b = 2.d0 !< [u]

real*8 :: x(2), s, R, Z, R_s, R_t, Z_s, Z_t, c
real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
integer :: i_elm

integer :: i, j, k, l, n_steps, n_lost, n_lost_all, i_elm_old, ifail, n_stream, seed, i_rng, ierr
real*8 :: target_time, t, n_e, T_e, grad_T_e(3), ran(6), q(3), coulomb_log, kTb, n_b, v_b(3,n_coll), ran2(6,n_coll)
!$ real*8 :: w0, w1
real*8 :: E(3), B(3), rz_old(2), st_old(2), psi, U, mmm(3)

real*8 :: proj_Lz_equil2

type(event) :: fieldreader
type(pcg32_rng), dimension(:), allocatable :: rng

type(type_edge_domain), allocatable, dimension(:) :: W_edge_domains, W_edge_domains2
type(edge_elements) :: W_edge, W_edge2

logical, parameter :: restart = .true.

! Start up MPI, jorek
call sim%initialize(num_groups=2)


! Prepare the coronal equilibrium
  

!========================INITIALISE PARTICLES AS FROM W_EXAMPLE==============================
if (.not. restart) then
  ! Set up the field reader
  fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
  !fieldreader = event(read_jorek_fields_interp_hermite_birkhoff(basename='jorek', i=1004))
  call with(sim, fieldreader)

  ! Set up particles
  sim%groups(1)%Z    = 74
  sim%groups(1)%mass = atomic_weights(74) !< atomic mass units
  sim%groups(1)%ad = read_adf11(sim%my_id,'50_w')
  sim%groups(1)%cor = coronal(sim%groups(1)%ad)
  sim%groups(2)%Z    = 74
  sim%groups(2)%mass = atomic_weights(74) !< atomic mass units
  sim%groups(2)%ad = sim%groups(1)%ad
  sim%groups(2)%cor = sim%groups(1)%cor


  do i=1,2
    allocate(particle_kinetic_leapfrog::sim%groups(i)%particles(n_particles/sim%n_cpu))

    !sim%groups(i)%particles%i_elm = 0 
    ! Uniformly initialise 10% of the particles
    !call initialise_particles_H_mu_psi(sim%groups(i)%particles(1:n_particles/(sim%n_cpu)), &
        !sim%fields, pcg32_rng(), &
        !sim%groups(i)%mass, &
        !uniform_space=.true., cor=sim%groups(i)%cor)

    select type (p => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      p(:)%i_elm = 0
      p(:)%weight = 5.d11 !< 1e5 particles uniformly over ~100 m^3 volume
      p(:)%t_birth = sim%time
      ! at 6.5e19 central density -> approx uniform concentration of 2.5e14 /m^3 at
      ! relative concentration of 1e-5
      ! -> 1e16 divided by 1e5 = 1e11
      call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(i)%mass, &
          sim%fields, sim%time, timesteps(i))
    end select
  end do
else
  call read_simulation_hdf5(sim, 'part_restart.h5')
  ! Set up the field reader
  !fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
  !fieldreader = event(read_jorek_fields_interp_hermite_birkhoff(basename='jorek', i=last_file_before_time(sim%time)))
  ! read this specific file and reset time
  sim%time = 0.d0 ! to cause the restart file reader to overwrite it
  fieldreader = event(read_jorek_fields_interp_hermite_birkhoff(basename='jorek', i=1002))
  call with(sim, fieldreader)
end if


! number of particles to sputter per species (should be renormalized to yield)
call find_edge_domains(sim%fields%node_list, sim%fields%element_list, W_edge_domains, sides=[.false., .false., .true., .false.]) ! inner div
call W_edge%prepare(sim%fields%node_list, sim%fields%element_list, W_edge_domains, nsub=8, nsub_toroidal=1)
! target group, number of particles per mpi task, densities, Zs, basename
W_sputter_source_in = particle_sputter(W_edge, 1, nint(real(n_particles)*sputter_dt/sputter_fill_time/real(sim%n_cpu)), [2.0*0.4635,0.04,0.03,0.003], [-2,2,7,18], basename='W_inner')
W_sputter_source_in%n_save = 1
W_sputter_source_in%use_thompson = .true.
W_sputter_source_in%use_Yn_func = .false.
!W_sputter_source_in%sputtered_particle_weight_threshold = 1d99 ! no cascade self-sputtering

call find_edge_domains(sim%fields%node_list, sim%fields%element_list, W_edge_domains2, sides=[.true., .false., .false., .false.]) ! outer div
call W_edge2%prepare(sim%fields%node_list, sim%fields%element_list, W_edge_domains2, nsub=8, nsub_toroidal=1)
W_sputter_source_out = particle_sputter(W_edge2, 2, nint(real(n_particles)*sputter_dt/sputter_fill_time/real(sim%n_cpu)), [2.0*0.4635,0.04,0.03,0.003], [-2,2,7,18], basename='W_outer')
W_sputter_source_out%n_save = 1
W_sputter_source_out%use_thompson = .true.
W_sputter_source_out%use_Yn_func = .false.
!W_sputter_source_out%sputtered_particle_weight_threshold = 1d99 ! no cascade self-sputtering

events = [fieldreader, &
          event(W_sputter_source_in, step=sputter_dt), &
          event(W_sputter_source_out, step=sputter_dt), & ! write diagnostics
          event(count_action(), step=1d-6), &
          !event(write_particle_diagnostics(filename='diag.h5', append=restart), step=1d-5), &
          event(write_action(),   step=1d-5), &
          event(projection(sim%fields%node_list, sim%fields%element_list, &
            filter=4d-5, filter_hyper=1d-9, f=[&
            proj_f(proj_one, 1), proj_f(proj_one, 2), proj_f(proj_Lz, 1), &
            proj_f(proj_Lz, 2), proj_f(proj_lz_equil, 1), proj_f(proj_Lz_equil, 2)],&
            to_h5=.true., basename='proj'),  step=1d-6), &
          event(stop_action(), start=1d-1)  &
        ]

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
    n_steps = floor((target_time - sim%time - step_rest_time(i))/timesteps(i))
    step_rest_time(i) = mod(target_time - sim%time - step_rest_time(i), timesteps(i))
    n_lost = 0

    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)

#ifdef __GFORTRAN__
      !$omp parallel do default(shared) & ! workaround for Error: ‘__vtab_mod_pcg32_rng_Pcg32_rng’ not specified in enclosing ‘parallel’
#else
      !$omp parallel do default(none) &
#endif
      !$omp shared(sim, n_steps, timesteps, i, rng, central_density) &
      !$omp private(i_rng, j, k, t, E, B, psi, U, rz_old, st_old, i_elm_old, n_e, &
      !$omp         T_e, grad_T_e, ran, ran2, kTb, n_b, q, coulomb_log, l, P, P_s, P_t, &
      !$omp         P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t, v_b, ifail) &
      !$omp schedule(static,10) ! Using dynamic with lower values can lead to inefficiences
      do j=1,size(particles,1)
        i_rng = 1
        !$ i_rng = omp_get_thread_num()+1

        do k=1,n_steps
          if (particles(j)%i_elm .le. 0) exit
          t = sim%time + k*timesteps(i)
          i_elm_old = particles(j)%i_elm
          call sim%fields%calc_EBpsiU(t, particles(j)%i_elm, &
              particles(j)%st, particles(j)%x(3), E, B, psi, U)
          rz_old    = particles(j)%x(1:2)
          st_old    = particles(j)%st
          i_elm_old = particles(j)%i_elm

          ! Update charge based on ionisation coefficients
          call sim%fields%calc_NeTe(t, i_elm_old, particles(j)%st, particles(j)%x(3), n_e, T_e, grad_T_e)
          call rng(i_rng)%next(ran)
          particles(j)%q = int(new_charge(int(particles(j)%q,4), sim%groups(i)%ad, log10(n_e), log10(T_e), timesteps(i), ran(1:2)),1)

          if (particles(j)%q .gt. 0) then
            ! Calculate collisions
            kTb = T_e*K_BOLTZ/EL_CHG ! assume T_e == T_i
            n_b = n_e
            q = q_homma2013(kTb, grad_T_e*EL_CHG/K_BOLTZ, B, n_b, m_b, q_b)
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
      !$omp end parallel do
    end select
  end do ! groups
  sim%time = target_time
  !$ w1 = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w1-w0)
  !$ if (sim%my_id .eq. 0) write(*,"(f10.7,A,3f9.4,A)") sim%time, " Particle stepping complete in ", mmm, "s"
  call with(sim, events, at=sim%time)
end do
call sim%finalize
end program W_sputtering_example
