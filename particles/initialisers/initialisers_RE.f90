!> module containing initialization and pdf (probability distribution function) generation
!> subroutines specifically useful for runaway electrons
module initialisers_RE
  use mod_particle_types
  use mod_particle_sim
  use mod_rng
  use initialisers_base
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT, SPEED_OF_LIGHT, MASS_ELECTRON, TWOPI, MU_ZERO
  use mod_pusher_tools, only: get_orthonormals
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_model_settings
  use phys_module, only: CENTRAL_DENSITY, CENTRAL_MASS
  use equil_info
  use mod_rej_f
  implicit none

  contains
subroutine initialise_re_gaussian(sim, group_num, rng, space_pdf, energy, pitch, std_energy)
  use phys_module, only: tstep_particles
  use mod_kinetic_relativistic
  use mod_sampling, only: boxmueller_transform

  type(particle_sim),                   intent(inout) :: sim
  integer,                              intent(in)    :: group_num
  class(type_rng),                      intent(in)    :: rng
  type(spatial_pdf),    optional,       intent(in)    :: space_pdf
  real*8,                               intent(in)    :: energy, pitch ! Kinetic energy in units of eV and pitch
  real*8,               optional,       intent(in)    :: std_energy
  real*8,               allocatable                   :: p_tot(:), p_par(:), p_perp(:)
  integer                                             :: j
  real(kind=8)                                        :: psi, U, gyro_angle
  real(kind=8),         dimension(3)                  :: E, B, B_cart, B_norm
  real*8                                              :: e1(3), e2(3) 
  integer                                             :: num_part
  real*8,               allocatable                   :: ran_uniform(:), ran_gaussian(:)

  call initialise_particles(sim%groups(group_num)%particles, &
    sim%fields%node_list, sim%fields%element_list, rng, space_pdf)

  num_part = size(sim%groups(group_num)%particles,1) 

  allocate(p_tot(num_part))
  allocate(p_par(num_part))
  allocate(p_perp(num_part))

  ! Generate gassian distributed energy
  if (present(std_energy)) then

    allocate(ran_uniform(num_part + mod(num_part,2)))
    allocate(ran_gaussian(num_part + mod(num_part,2)))

    call random_number(ran_uniform)
    ran_gaussian = boxmueller_transform(ran_uniform)

    p_tot               = sqrt(((energy + std_energy*ran_gaussian(:num_part))*EL_CHG/SPEED_OF_LIGHT + MASS_ELECTRON*SPEED_OF_LIGHT)**2 - (MASS_ELECTRON*SPEED_OF_LIGHT)**2)/ATOMIC_MASS_UNIT ! [AMU*m/s]

    deallocate(ran_uniform)
    deallocate(ran_gaussian)

  else

    p_tot               = sqrt((energy*EL_CHG/SPEED_OF_LIGHT + MASS_ELECTRON*SPEED_OF_LIGHT)**2 - (MASS_ELECTRON*SPEED_OF_LIGHT)**2)/ATOMIC_MASS_UNIT ! [AMU*m/s]

  end if

  p_par               = pitch * p_tot
  p_perp              = sqrt(p_tot**2 - p_par**2)

  ! Set particle momentum
  select type (particles => sim%groups(group_num)%particles)
  type is (particle_kinetic_relativistic)
    !$omp parallel do default(none) &
    !$omp private(E, B, psi, U, B_cart, B_norm, e1, e2, gyro_angle, j) &
    !$omp shared (sim, tstep_particles, p_par, p_perp, num_part)
    do j=1,num_part

      ! Extract magnetic field and convert to cartesian coordinates
      call sim%fields%calc_EBpsiU(sim%time, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
      B_cart = vector_cylindrical_to_cartesian(particles(j)%x(3),B)

      B_norm = B_cart/norm2(B_cart)

      ! Generate perpendicular component based on sampled gyro angle
      call get_orthonormals(B_norm, e1, e2)
      call random_number(gyro_angle)
      gyro_angle = gyro_angle * TWOPI

      particles(j)%p = p_par(j) * B_norm + p_perp(j)*(e1*cos(gyro_angle) + e2*sin(gyro_angle))

    end do
    !$omp end parallel do 
  end select

  deallocate(p_tot)
  deallocate(p_par)
  deallocate(p_perp)

end subroutine initialise_re_gaussian

end module initialisers_RE
