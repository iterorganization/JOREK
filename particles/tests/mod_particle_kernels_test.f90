#ifdef CUDA_KERNELS
module mod_particle_kernels_test
  use fruit
  use cudafor
  use mpi
  use mod_particle_kernels
  use mod_particle_types, only: particle_kinetic_leapfrog
  use mod_particle_sim, only: particle_sim
  implicit none
  
  private
  public :: run_fruit_particle_kernels

  type(particle_sim) :: sim

contains

  subroutine run_fruit_particle_kernels()
      write(*,'(/A)') "  ... setting-up: particle kernels tests"
      call setup
      write(*,'(/A)') "  ... running: particle kernels tests"
      call test_copy_device_data
      call test_particle_kinetic_leapfrog_loop
      write(*,'(/A)') "  ... tearing-down: particle kernels tests"
      call teardown
  end subroutine run_fruit_particle_kernels

  !> Set-up and tear-down -------------------------
  !> set-up particle types test features
  subroutine setup()

    !> Initialise a particle sim, we shouldn't need this but it is an easy way to set
    !> up the fields that we do need. We will ignore the particle group within the sim.
    call sim%initialize(num_groups=1)
    
    !> Set up the fields from a restart file. This would ideally be done in situ by calling a
    !> fields setup written for fields unti tests and not use a restart file.
    !> This approach will suffice for now but should remain under review
    call init_fields(sim)

    !> Set up particles using the H_mu_psi initialisation - this is not ideal but sufficient
    call init_particles(sim)

  end subroutine setup

  subroutine teardown()
  end subroutine teardown


  !> Reads fields from an hdf5 restart file jorek_restart.h5
  subroutine init_fields(sim_in)
    
    use mod_event
    use mod_fields_linear
    use equil_info,               only: update_equil_state
    use phys_module,              only: xcase, xpoint
    use mod_boundary,             only: boundary_from_grid
    !    use data_structure, only: type_bnd_element_list, type_bnd_node_list 
    use nodes_elements,           only: bnd_node_list, bnd_elm_list
    
    !> variables
    type(event)                       :: fieldreader
    type(particle_sim),intent(inout)  :: sim_in
    
    write(*,*) "Initialising fields from hdf5 restart file"
    
    ! Set up a field reader
    fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
    call with(sim_in, fieldreader)
    
    if (sim_in%my_id .eq. 0) call boundary_from_grid(sim_in%fields%node_list, sim_in%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
    
    call broadcast_boundary(sim_in%my_id, bnd_elm_list, bnd_node_list)
    
    call update_equil_state(sim_in%my_id, sim_in%fields%node_list, sim_in%fields%element_list, bnd_elm_list, xpoint, xcase)
    
  end subroutine init_fields


  !> Initialise particles on the existing geometry using initialise_particles_H_mu_psi
  subroutine init_particles(sim_in)

    use mod_atomic_elements,      only: atomic_weights
    use mod_initialise_particles, only: initialise_particles_H_mu_psi, adjust_particle_weights
    use mod_boris,                only: boris_all_initial_half_step_backwards_RZPhi
    use phys_module,              only: n_particles
    use mod_particle_types,       only: particle_kinetic_leapfrog
    use mod_pcg32_rng

    type(particle_sim),intent(inout)  :: sim_in
    type(pcg32_rng), dimension(:), allocatable     :: rng

    real*8                                         :: rho_part, timesteps
    integer                                        :: n_particles_local

    write(*,*) "Initialising particles using random distribution"

    ! Set up particles
    sim_in%groups(1)%Z    = 1
    sim_in%groups(1)%mass = atomic_weights(-2) !< atomic mass units

    rho_part           = 1.195d19
    timesteps          = 1e-10

    n_particles_local = int(n_particles/sim_in%n_cpu) 
    allocate(particle_kinetic_leapfrog::sim_in%groups(1)%particles(n_particles_local))

    select type (p => sim_in%groups(1)%particles)
    type is (particle_kinetic_leapfrog)

       call initialise_particles_H_mu_psi(p, sim_in%fields, pcg32_rng(),sim_in%groups(1)%mass, &
            uniform_space=.true., uniform_space_rej_f=f_toroidal_flux, &
            uniform_space_rej_vars=[1], charge = 1, T_maxwell = 4d5)

       call adjust_particle_weights(sim_in%groups(1)%particles, rho_part)
       if (sim_in%my_id .eq. 0) write(*,*) "Particle density was adjusted to:", rho_part, sim_in%groups(1)%particles(1:10)%weight

       call boris_all_initial_half_step_backwards_RZPhi(p, sim_in%groups(1)%mass, sim_in%fields, sim_in%time, timesteps)

       write(*,*) "Initialised particles",size(sim_in%groups(1)%particles,1)

    end select

  end subroutine init_particles

  subroutine test_copy_device_data

    type(particle_kinetic_leapfrog), managed, dimension(:), allocatable  :: particles
    type(fields_linear_device), managed, allocatable                      :: fields

    integer    :: n_particles

    call copy_device_data( sim , particles , fields )

  end subroutine test_copy_device_data
  
  
  subroutine test_particle_kinetic_leapfrog_loop

    call particle_kinetic_leapfrog_loop( sim )

    
  end subroutine test_particle_kinetic_leapfrog_loop    


  !> required for particle initialisation
  pure function f_toroidal_flux(n, P, grad_P) result(f)

    use equil_info,               only: ES

    integer, intent(in) :: n
    real*8, intent(in) :: P(n), grad_P(3,n)
    real*8 :: s, psi_norm, coeff(0:3)
    real*4 :: f

    ! central densiy should be 1.44131x10^17

    coeff(0)=0.49123
    coeff(1)=0.298228
    coeff(2)=0.198739
    coeff(3)=0.521298

    psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

    s = 0.957 * psi_norm + 0.043 * psi_norm**2 

    f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

  end function f_toroidal_flux
 
end module mod_particle_kernels_test
#endif
