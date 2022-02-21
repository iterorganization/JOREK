!> Module for the initialization of multiple relativistic
!> particles in the phase space (3D spatial 2D/3D 
!> velocity coordinates. Acceptance-rejection method
!> is used for the Monte-Carlo initialisation
module mod_initialise_relativistic_particles
implicit none
private
#ifdef UNIT_TESTS
public :: sample_position_uniformly_cylinder
public :: sample_position_uniformly_psi_theta
public :: sampling_cartesian_p_kinetic_relativistic
public :: sampling_cartesian_p_gc_relativistic
public :: sampling_uniform_ppitchgyro_kinetic_relativistic
public :: sampling_uniform_ppitchgyro_gc_relativistic
#endif

!> Variables and datatypes --------------------------------
!> Interfaces ---------------------------------------------
interface
  !> interface of the find RZ procedure
  subroutine find_RZ(node_list,element_list,R_find,Z_find,&
  R_out,Z_out,ielm_out,s_out,t_out,ifail)
    use data_structure
    use mod_element_rtree, only: elements_containing_point
    implicit none
    type(type_node_list),intent(in)    :: node_list
    type(type_element_list),intent(in) :: element_list
    real*8,intent(in)                  :: R_find,Z_find
    real*8,intent(out)                 :: R_out,Z_out,s_out,t_out
    integer,intent(inout)              :: ielm_out
    integer,intent(out)                :: ifail
  end subroutine find_RZ

  !> interface of the find theta psi procedure
  subroutine find_theta_psi(node_list,element_list,psi_minmax,&
  theta,psi,phi,R_axis,Z_axis,i_elm,s,t,R,Z)
    use data_structure
    implicit none
    !> inputs:
    type(type_node_list),intent(in)    :: node_list
    type(type_element_list),intent(in) :: element_list
    real*8,dimension(element_list%n_elements,2),intent(in) :: psi_minmax
    real*8,intent(in)                  :: theta,psi,phi
    real*8,intent(in)                  :: R_axis,Z_axis
    integer,intent(in)                 :: i_elm
    !> outputs:
    real*8,intent(out)                 :: s,t,R,Z 
  end subroutine find_theta_psi

  !> find minim and maximum value of the poloidal flux (psi)
  !> for each mesh element
  subroutine psi_minmax(node_list,element_list,i_elm,psimin,psimax)
    use data_structure
    implicit none
    type(type_node_list)    :: node_list
    type(type_element_list) :: element_list
    real*8                  :: psimin,psimax
    integer                 :: i_elm
  end subroutine psi_minmax
end interface
contains
!> Procedures ---------------------------------------------
!> sample particles in the physical space. Particles
!> are distributed uniformely in the physical space
!> using the cylindrical coordinates.
!> inputs:
!>   particle:     (particle_base) particle to be sampled
!>   node_list:    (type_node_list) jorek element nodes
!>   element_list: (type_element_list) jorek element list
!>   rand:         (real8)(3) random numbers
!>   R2_bounds:    (real8)(2) squared of the major radius 
!>                            bounding box
!>   Z_bound:      (real8)(2) vertical position bounding box
!>   phi_bound:    (real8)(3) toroidal angle bounding box
!> outputs:
!>   particle: (particle_base) sampled particle
subroutine sample_position_uniformly_cylinder(particle,&
node_list,element_list,rand,R2_bound,Z_bound,phi_bound)
  use data_structure,     only: type_node_list
  use data_structure,     only: type_element_list
  use mod_particle_types, only: particle_base
  use mod_sampling,       only: transform_uniform_cylindrical
  implicit none
  !> inputs-outputs:
  class(particle_base),intent(inout) :: particle
  !> inputs:
  type(type_node_list),intent(in)    :: node_list
  type(type_element_list),intent(in) :: element_list
  real*8,dimension(2),intent(in)     :: R2_bound,Z_bound,phi_bound
  real*8,dimension(3),intent(in)     :: rand
  !> variables:
  integer :: ifail
  !> compute the particle position in global and local coords
  call transform_uniform_cylindrical(rand,R2_bound,Z_bound,&
  phi_bound,particle%x)
  call find_RZ(node_list,element_list,particle%x(1),particle%x(2),&
  particle%x(1),particle%x(2),particle%i_elm,particle%st(1),&
  particle%st(2),ifail)
end subroutine sample_position_uniformly_cylinder

!> sample particles in physical space. Particles are 
!> distributed uniformly in the poloidal flux and 
!> poloidal angle coordinates which implies that 
!> the distribution is not uniformly sampled in in 
!> physical space
!> inputs:
!>   particle:           (particle_base) particle to be sampled
!>   node_list:          (type_node_list) jorek element nodes
!>   element_list:       (type_element_list) jorek element list
!>   rand:               (real8)(3) random numbers
!>   thetapsiphi_bound:  (real8)(3,2) sampling box boundaries
!>                       1) minimum theta,psi,phi boundary
!>                       2) maximum theta,psi,phi boundary
!>   psi_element_minmax: (real8)(n_elements,2) minimum and 
!>                       maximum value of psi per each element
!>   RZ_axis:            (real8)(2) 1)-R 2)-Z magnetic axis coords.
!> outputs: 
!>   particle:           (particle_base) sampled particle
subroutine sample_position_uniformly_psi_theta(&
particle,node_list,element_list,rand,&
thetapsiphi_bound,psi_element_minmax,RZ_axis)
  use data_structure,     only: type_node_list
  use data_structure,     only: type_element_list
  use mod_particle_types, only: particle_base
  implicit none
  !> inputs-outputs:
  class(particle_base),intent(inout) :: particle
  !> inputs:
  type(type_node_list),intent(in)    :: node_list
  type(type_element_list),intent(in) :: element_list
  real*8,dimension(2),intent(in)     :: RZ_axis
  real*8,dimension(3),intent(in)     :: rand
  real*8,dimension(3,2),intent(in)   :: thetapsiphi_bound
  real*8,dimension(element_list%n_elements,2),intent(in) :: psi_element_minmax
  !> compute particle position in global coordinates
  particle%x = thetapsiphi_bound(:,1)+(thetapsiphi_bound(:,2)-&
  thetapsiphi_bound(:,1))*rand
  call find_theta_psi(node_list,element_list,psi_element_minmax,&
  particle%x(1),particle%x(2),particle%x(3),RZ_axis(1),RZ_axis(2),&
  particle%i_elm,particle%st(1),particle%st(2),particle%x(1),particle%x(2))
end subroutine sample_position_uniformly_psi_theta

!> uniform sampling of the cartesian momentum
!> for kinetic relativistic particles
!> inputs:
!>   particle:        (particle_kinetic_relativistic) particle to be sampled
!>   rand:            (real8)(2) random numbers
!>   pxpypz_interval: (real8)(3,2) momenta interval (in units of AMU)
!>                    1- px,py,pz momentum lower intervals
!>                    2- px,py,pz momentum upper intervals
!> outputs:
!>   particle: (particle_base) particle with sample momenta
subroutine sampling_cartesian_p_kinetic_relativistic(particle,&
rand,pxpypz_interval)
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  !> inputs-outputs:
  type(particle_kinetic_relativistic),intent(inout) :: particle
  !> inputs:
  real*8,dimension(3),intent(in)   :: rand
  real*8,dimension(3,2),intent(in) :: pxpypz_interval
  !> sample the cartesian momentum
  !> WARNING: the momentum must be in unit of AMU
  particle%p = pxpypz_interval(:,1)+(pxpypz_interval(:,2)-pxpypz_interval(:,1))*rand
end subroutine sampling_cartesian_p_kinetic_relativistic

!> uniform sampling of the cartesiam momentum 
!> for gc relativistic particles
!> inputs:
!>   gc:              (particle_gc_relativistic) relativistic gc to be sampled
!>   fields:          (fields_base) MHD fields type object
!>   rand:            (real8)(3) random numbers
!>   time:            (real8) time at which perform the interpolation
!>   mass:            (real8) particle mass in AMU
!>   pxpypz_interval: (real8)(3,2) momenta interval (in units of AMU)
!>                    1- px,py,pz momentum lower intervals
!>                    2- px,py,pz momentum upper intervals
!> outputs:
!>   gc_inout: (particle_gc_relativistic) relativistic gc with sampled momentum
subroutine sampling_cartesian_p_gc_relativistic(gc,fields,&
rand,time,mass,pxpypz_interval)
  use mod_particle_types,       only: particle_base
  use mod_particle_types,       only: particle_gc_relativistic
  use mod_fields,               only: fields_base
  use mod_kinetic_relativistic, only: momentum_relativistic_kinetic_to_relativistic_gc
  implicit none
  !> inputs-outputs:
  type(particle_gc_relativistic),intent(inout) :: gc
  class(fields_base),intent(inout) :: fields
  !> inputs:
  real*8,intent(in)                :: time,mass
  real*8,dimension(3),intent(in)   :: rand
  real*8,dimension(3,2),intent(in) :: pxpypz_interval
  !> variables:
  real*8 :: psi,U,B_norm 
  real*8,dimension(3) :: momentum,B,E
  !> sample the relativistic gc
  !> WARNING: the momentum must be in units of AMU
  momentum = pxpypz_interval(:,1)+(pxpypz_interval(:,2)-pxpypz_interval(:,1))*rand
  call fields%calc_EBpsiU(time,gc%i_elm,gc%st,gc%x(3),&
  E,B,psi,U)
  B_norm = sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3)); B = B/B_norm;
  gc%p = momentum_relativistic_kinetic_to_relativistic_gc(mass,&
  momentum,gc%x(3),B_norm,B)
end subroutine sampling_cartesian_p_gc_relativistic

!> uniform sampling of the energy, pitch angle
!> and gyro angle for kinetic relativistic particles
!> inputs:
!>   particle:       (particle_kinetic_relativistic) particle to be sampled
!>   fields:         (fields_base) MHD fields type object
!>   rand:           (real8)(3) random numbers
!>   time:           (real8) time of the initialisation
!>   R_cube_int:     (real8)(2) cube of the major radius sampling interval
!>   cos_theta_int:  (real8)(2) cosinus of the pitch angle sampling interval
!>   gyro_int:       (real8)(2) gyro angle sampling interval
!> outputs:
!>   particle: (particle_kinetic_relativistic) particle with sampled momentum
subroutine sampling_uniform_ppitchgyro_kinetic_relativistic(particle,&
fields,rand,time,R_cube_int,cos_theta_int,gyro_int)
  use mod_sampling,             only: sample_uniform_sphere_corona_rcosphi
  use mod_particle_types,       only: particle_kinetic_relativistic
  use mod_fields,               only: fields_base
  use mod_kinetic_relativistic, only: kinetic_relativistic_momentum_spherical_to_cart
  implicit none
  !> inputs-outputs:
  type(particle_kinetic_relativistic),intent(inout) :: particle
  class(fields_base),intent(inout) :: fields
  !> inputs:
  real*8,intent(in)               :: time
  real*8,dimension(2),intent(in)  :: R_cube_int,cos_theta_int,gyro_int
  real*8,dimension(3),intent(in)  :: rand
  !> variables:
  real*8              :: psi,U
  real*8,dimension(3) :: B,E,rpitchgyro
  !> compute the magnetic field direction
  call fields%calc_EBpsiU(time,particle%i_elm,particle%st,&
  particle%x(3),E,B,psi,U)
  B = B/sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))
  !> sample the momentum in spherical coordinates
  rpitchgyro = sample_uniform_sphere_corona_rcosphi(R_cube_int,&
  cos_theta_int,gyro_int,rand)
  !> compute the particle kinetic momentum
  particle%p = kinetic_relativistic_momentum_spherical_to_cart(&
  particle%x(3),particle%p,B)
end subroutine sampling_uniform_ppitchgyro_kinetic_relativistic

!> uniform sampling of the energy, pitch angle
!> and gyro angle for gc relativistic particles
!> inputs:
!>   particle:       (particle_gc_relativistic) particle to be sampled
!>   fields:         (fields_base) MHD fields type object
!>   rand:           (real8)(3) random numbers
!>   mass:           (real8) particle mass in AMU
!>   time:           (real8) time of the initialisation
!>   R_cube_int:     (real8)(2) cube of the major radius sampling interval
!>   cos_theta_int:  (real8)(2) cosinus of the pitch angle sampling interval
!>   gyro_int:       (real8)(2) gyro angle sampling interval
!> outputs:
!>   particle: (particle_gc_relativistic) particle with sampled momentum
subroutine sampling_uniform_ppitchgyro_gc_relativistic(gc,fields,&
rand,mass,time,R_cube_int,cos_theta_int,gyro_int)
  use mod_sampling,             only: sample_uniform_sphere_corona_rcosphi
  use mod_particle_types,       only: particle_gc_relativistic
  use mod_fields,               only: fields_base
  use mod_kinetic_relativistic, only: kinetic_relativistic_momentum_spherical_to_cart
  use mod_kinetic_relativistic, only: momentum_relativistic_kinetic_to_relativistic_gc
  implicit none
  !> inputs-outputs:
  type(particle_gc_relativistic),intent(inout) :: gc
  class(fields_base),intent(inout) :: fields
  !> inputs:
  real*8,intent(in)               :: time,mass
  real*8,dimension(2),intent(in)  :: R_cube_int,cos_theta_int,gyro_int
  real*8,dimension(3),intent(in)  :: rand
  !> variables
  real*8              :: B_norm,psi,U
  real*8,dimension(3) :: p_kin,B,E
  !> compute magnetic field direction
  call fields%calc_EBpsiU(time,gc%i_elm,gc%st,gc%x(3),&
  E,B,psi,U)
  B_norm = sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3)); B = B/B_norm;
  !> sample the momentum in spherical coordinates
  p_kin = sample_uniform_sphere_corona_rcosphi(R_cube_int,&
  cos_theta_int,gyro_int,rand)
  !> compute the particle kinetic momentum
  p_kin = kinetic_relativistic_momentum_spherical_to_cart(gc%x(3),p_kin,B)
  !> transform in guiding center coords
  gc%p = momentum_relativistic_kinetic_to_relativistic_gc(mass,&
  p_kin,gc%x(3),B_norm,B)
end subroutine sampling_uniform_ppitchgyro_gc_relativistic

!>---------------------------------------------------------
end module mod_initialise_relativistic_particles
