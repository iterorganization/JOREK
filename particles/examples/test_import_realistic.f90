program test_import_realistic

  ! Import everything
  use particle_tracer
  use mod_particle_diagnostics
  use mpi
  use mod_atomic_elements
  use mod_particle_io
  use mod_event
  use mod_project_particles
  use mod_particle_loop
  use mod_jorek_timestepping
  use mod_random_seed
  use mod_interp, only: mode_moivre, interp_RZ, interp_0
  use mod_basisfunctions
  use nodes_elements
  use phys_module, only: tstep, restart, t_start, restart_particles
  use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint
  use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
  use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
  use phys_module, only: n_mode_families

  use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG

  use mod_edge_domain
  use mod_edge_elements
  use data_structure, only: type_bnd_element_list, type_bnd_node_list
  use equil_info
  use mod_boundary, only: boundary_from_grid

  ! Import realistic distribution function
  use mod_import_experimental_dist
  !use mod_phase_space_project ! Used for looking at 4D distribution function, but not in develop. Commented out.
!$ use omp_lib
  implicit none

  integer                                               :: n_particles_local,ierr
  type(event)                                           :: fieldreader,partwriter
  !real*8, allocatable                                  :: phase_proj(:)
  !type(phase_space_projection)                         :: test_phase
  type(particle_kinetic_leapfrog)                       :: particle_tmp
  type(projection), target                              :: jorek_feedback
  real*8                                                :: timesteps,E(3),B(3),psi,U,Energy,t, rho_part
  !real*8, dimension(:),allocatable                     :: val_tmp
  !integer, dimension(:),allocatable                    :: index_phase_tmp
  integer                                               :: i_phase,j, i_elm,i_tor,l,m, index_lm
  real*8,allocatable                                    :: feedback_rhs(:,:,:,:,:)
  real*8                                                :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4),v
  
  write(*,*) "Starting test of importing a realistic distribution function..."
  call sim%initialize(num_groups=1)
  rho_part    = 1.195d19*2.d0
  n_particles_local = int(n_particles/sim%n_cpu)
  timesteps         = tstep_particles
  ! Read JOREK fields
  fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
  call with(sim, fieldreader)
  if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
  call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)
  call update_equil_state(sim%my_id,sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)


  ! Set up particles
  sim%groups(1)%Z    = 1
  sim%groups(1)%mass = atomic_weights(-2) !< atomic mass units
  allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(n_particles_local))

  call import_particles(sim%groups(1)%particles,sim%fields,"experimental_dist.h5",pcg32_rng(),sim%groups(1)%mass,4,0.999d0)

  select type (p => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)

      call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(1)%mass, sim%fields, sim%time, timesteps)

  end select

  call adjust_particle_weights(sim%groups(1)%particles, rho_part)

  partwriter = event(write_action(filename='part_restart.h5'))
  call with(sim, partwriter)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !test_phase = new_phase_space_projection(ndim=3,res=[200,200,100],start=[0d0,-1d0,0.d0],end=[100d0,1.1d0,1.d0],f_proj=proj_f(proj_one, group = 1),&
  !                                       f_grids= [proj_f(proj_one,group = 1 ),proj_f(proj_one,group = 1 )],bandwidths=[5d0,0.3d0,0.2d0])

  jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                         filter    = filter_perp, filter_hyper    = filter_hyper, filter_parallel    = filter_par, &
                                         filter_n0 = filter_perp, filter_hyper_n0 = filter_hyper, filter_parallel_n0 = filter_par_n0, &
                                         calc_integrals=.false., to_vtk=.true., to_h5 = .false., basename='projections')
  allocate(jorek_feedback%rhs(n_degrees, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 1))

  jorek_feedback%rhs = 0.d0
  allocate(feedback_rhs,source=jorek_feedback%rhs)                                         

  !allocate(phase_proj,source=test_phase%values)
  !allocate(val_tmp(test_phase%totsupport))
  !allocate(index_phase_tmp(test_phase%totsupport))
  !phase_proj(:) = 0.d0
  t=0.d0
  select type(particles=> sim%groups(1)%particles)
    type is(particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
      !$omp parallel do default(shared) & 
#else
      !$omp parallel do default(none)&
      !$omp shared(sim,particles,t)&
#endif
      !$omp private(j,particle_tmp,E,B,psi,U,Energy,i_elm,i_tor,l,m,index_lm,v,HH,HZ,HH_s,HH_t)&
      !$omp reduction(+:feedback_rhs)
      do j=1,size(particles)
        call copy_particle_kinetic_leapfrog(particles(j),particle_tmp)
        call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
        Energy =0.5d0*sim%groups(1)%mass*MASS_PROTON*dot_product(particle_tmp%v,particle_tmp%v)/EL_CHG/1000d0
        !call calc_index_shaped_part_x(test_phase,index_phase_tmp, val_tmp,[Energy,dot_product(particle_tmp%v,B)/norm2(B)/norm2(particle_tmp%v),sqrt((psi+0.16192)/(0.16192))])
        !do i_phase=1, test_phase%totsupport
        !  if(index_phase_tmp(i_phase) > 0)then
        !    phase_proj(index_phase_tmp(i_phase))=phase_proj(index_phase_tmp(i_phase))+1.d0*particle_tmp%weight*val_tmp(i_phase)!*E_diff
        !
        !  endif
        !enddo

        i_elm = particle_tmp%i_elm

        if (i_elm .gt. 0) then

          call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
          call mode_moivre(particle_tmp%x(3), HZ)

          do l=1,n_vertex_max
            do m=1,n_degrees

              index_lm = (l-1)*n_degrees + m

              v = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m)

              do i_tor=1,n_tor
                feedback_rhs(m,l,i_elm,i_tor,1) = feedback_rhs(m,l,i_elm,i_tor,1) &

                                                  + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                  * (1.d0/3.d0) * (particle_tmp%v(1)**2 + particle_tmp%v(2)**2 + particle_tmp%v(3)**2) * mu_zero
              enddo !< n_tor

            enddo   !< order
          enddo     !< vertex
        endif
      enddo

      !$omp end parallel do
      jorek_feedback%rhs = feedback_rhs
      !test_phase%values = phase_proj
      !call output_phase_project(test_phase,1)
      call with(sim,jorek_feedback)
  end select

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_ABORT(MPI_COMM_WORLD,10,ierr)
  
  
end program test_import_realistic
