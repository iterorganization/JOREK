program test_generalised_initialisation_gc_H_mu_psi
!> the program is for comparing the generic initialisation method 
!> contained in mod_initialise_particles and the case specific 
!> initialize_particles_H_mu_psi implemented in mod_initialise_particles
use constants, only: TWOPI,PI,ATOMIC_MASS_UNIT,EL_CHG
use particle_tracer
use mod_interp
implicit none
type(pcg32_rng) :: rng_pcg32
type(event)     :: field_reader,particle_writer
integer         :: n_variables,n_particles,n_int_pdf_param,n_real_pdf_param,ifail
integer         :: n_int_weight_param,n_real_weight_param
integer         :: n_int_gdf_param,n_real_gdf_param
integer         :: n_int_pdf_to_part_coord_param,n_real_pdf_to_part_coord_param
integer         :: include_vpar,particle_type,t1,t2,t3,t4,dummy_int
integer,dimension(:),allocatable :: int_pdf_param,int_weight_param,int_gdf_param
integer,dimension(:),allocatable :: int_pdf_to_part_coord_param
real*8          :: start_time,mass,pdf_upper_bound,gdf_upper_bound
real*8          :: dummy_real8_1,dummy_real8_2,dummy_real8_3
real*8,dimension(2)                 :: Poloidalbound,Phibound 
real*8,dimension(2)                 :: Ekinbound,Vbound,Pitchbound,Gyrobound,Chargebound
real*8,dimension(:),allocatable     :: psi_minmax_loc,real_pdf_to_part_coord_param
real*8,dimension(:),allocatable     :: real_pdf_param,real_weight_param,real_gdf_param
real*8,dimension(7,2)               :: phase_space_bounds
character(len=125)                  :: generic_particle_filename,original_h_mu_psi_filename
character(len=:),allocatable        :: jorek_filename
procedure(real_f),pointer           :: pdf_to_use         => NULL()
procedure(real_f),pointer           :: weight_to_use      => NULL()
procedure(real_f),pointer           :: gdf_to_use         => NULL()
procedure(real_arr_inout_s),pointer :: gdf_sampler_to_use => NULL()
procedure(part_inout_s),pointer     :: pdf_to_part_coord  => NULL()

! MPI and group initialisation -------------------------------------------------------------------
call sim%initialize(num_groups=1)

! Define inputs -------------- -------------------------------------------------------------------
generic_particle_filename  = 'part_restart_H_mu_psi_generic.h5'
original_h_mu_psi_filename = 'part_restart_H_mu_psi_original.h5'
particle_type             = 1
n_variables               = 7
n_particles               = 100000000
include_vpar              = 0
start_time                = 0d0
mass                      = 2.01410177811*ATOMIC_MASS_UNIT 
Poloidalbound             = [0d0,TWOPI]
Phibound                  = [0d0,TWOPI]
Ekinbound                 = [-9.99d9,9.99d9] !< in eV not required
Pitchbound                = [0d0,PI]
Gyrobound                 = [0d0,TWOPI]
Chargebound               = 1d0
allocate(character(25) :: jorek_filename)
jorek_filename            = 'jorek_equilibrium'

!> Initialisations -------------------------------------------------------------------------------
write(*,*) "Test H_mu_psi initialisation generic vs original, initialise parameters ... "
field_reader = event(read_jorek_fields_interp_linear(basename=trim(jorek_filename),i=-1))
call with(sim,field_reader)
!> initailise simulation parameters
sim%time = start_time
sim%groups(1)%mass = mass
!> Define the sample box
Vbound = sqrt(2.d0*Ekinbound*EL_CHG/mass)
allocate(psi_minmax_loc(2*sim%fields%element_list%n_elements))
call extract_element_psi_minmax(sim%fields,psi_minmax_loc)
phase_space_bounds(:,1) = [minval(psi_minmax_loc(1:sim%fields%element_list%n_elements)),&
                          Poloidalbound(1),Phibound(1),Vbound(1),Pitchbound(1),Gyrobound(1),Chargebound(1)]
phase_space_bounds(:,2) = [maxval(psi_minmax_loc(sim%fields%element_list%n_elements+1:&
                          2*sim%fields%element_list%n_elements)),Poloidalbound(2),Phibound(2),&
                          Vbound(2),Pitchbound(2),Gyrobound(2),Chargebound(2)]
!> define the pdf inputs
pdf_to_use         => pdf_psi_H_mu
n_int_pdf_param    = 0
n_real_pdf_param   = 0
pdf_upper_bound    = sup_pdf_psi_H_mu(n_variables,&
phase_space_bounds(:,1),phase_space_bounds(:,2),&
n_real_pdf_param,real_pdf_param,&
n_int_pdf_param,int_pdf_param)
!> define the gdf inputs
gdf_to_use         => gdf_psi_H_mu
gdf_sampler_to_use => gdf_psi_H_mu_sampler
n_int_gdf_param    = 1
allocate(int_gdf_param(n_int_gdf_param))
int_gdf_param(1)   = include_vpar
n_real_gdf_param   = 2*sim%fields%element_list%n_elements+3
allocate(real_gdf_param(n_real_gdf_param))
real_gdf_param(1:2*sim%fields%element_list%n_elements) = psi_minmax_loc
deallocate(psi_minmax_loc)
call find_axis(sim%my_id,sim%fields%node_list,sim%fields%element_list,dummy_real8_1,&
real_gdf_param(2*sim%fields%element_list%n_elements+1),&
real_gdf_param(2*sim%fields%element_list%n_elements+2),dummy_int,dummy_real8_2,dummy_real8_3,ifail)
real_gdf_param(2*sim%fields%element_list%n_elements+3) = mass
gdf_upper_bound   = sup_gdf_psi_H_mu(n_variables,&
phase_space_bounds(:,1),phase_space_bounds(:,2),&
n_real_gdf_param,real_gdf_param,&
n_int_gdf_param,int_gdf_param)
!> define the weights inputs:
weight_to_use       => particle_weight_psi_H_mu
n_int_weight_param  = 0
n_real_weight_param = 0
!> define the sample to particle transformation inputs
pdf_to_part_coord  => copy_H_mu_psi_sample_to_p_gc
n_int_pdf_to_part_coord_param = 0
n_real_pdf_to_part_coord_param = 1
allocate(real_pdf_to_part_coord_param(n_real_pdf_to_part_coord_param))
real_pdf_to_part_coord_param(1) = mass

write(*,*) "Test H_mu_psi initialisation generic vs original, initialise parameters: completed"

write(*,*) "Test H_mu_psi initialisation generic vs original, initialise particles generic method ... "
call allocate_particle_list(particle_type,n_particles,sim%groups(1)%particles)
particle_writer = event(write_action(filename=trim(generic_particle_filename)))
call system_clock(t1)
call initialise_particles_in_phase_space(n_variables,sim%groups(1)%particles,sim%fields,rng_pcg32,&
pdf_to_use,weight_to_use,gdf_to_use,gdf_sampler_to_use,pdf_upper_bound,gdf_upper_bound,&
pdf_to_part_coord,sim%groups(1)%mass,start_time,phase_space_bounds,n_real_pdf_param,real_pdf_param,&
n_int_pdf_param,int_pdf_param,n_real_weight_param,real_weight_param,n_int_weight_param,&
int_weight_param,n_real_gdf_param,real_gdf_param,n_int_gdf_param,int_gdf_param,&
n_real_pdf_to_part_coord_param,real_pdf_to_part_coord_param,&
n_int_pdf_to_part_coord_param,int_pdf_to_part_coord_param)
call system_clock(t2)
call with(sim,particle_writer)
write(*,*) "Test H_mu_psi initialisation generic vs original, initialise particles generic method: completed!"

write(*,*) "Test H_mu_psi initialisation generic vs original, initialise particles originasl method ... "
call allocate_particle_list(particle_type,n_particles,sim%groups(1)%particles)
particle_writer = event(write_action(filename=trim(original_h_mu_psi_filename)))
call system_clock(t3)
call initialise_particles_H_mu_psi(sim%groups(1)%particles, sim%fields,rng_pcg32,mass,&
include_vpar=include_vpar.eq.1, charge=int(Chargebound(1)))
call system_clock(t4)
call with(sim,particle_writer)
write(*,*) "Test H_mu_psi initialisation generic vs original, initialise particles original method: completed!"
write(*,*) sim%my_id, 'System time H mu psi initialisation generic method (s): ',real(t2-t1,kind=8)/1d3
write(*,*) sim%my_id, 'System time H mu psi initialisation original method (s): ',real(t4-t3,kind=8)/1d3

! Clean-up ---------------------------------------------------------------------------------------
if(allocated(psi_minmax_loc))               deallocate(psi_minmax_loc)
if(allocated(int_pdf_param))                deallocate(int_pdf_param)
if(allocated(int_weight_param))             deallocate(int_weight_param)
if(allocated(int_gdf_param))                deallocate(int_gdf_param)
if(allocated(int_pdf_to_part_coord_param))  deallocate(int_pdf_to_part_coord_param)
if(allocated(real_pdf_to_part_coord_param)) deallocate(real_pdf_to_part_coord_param)
if(allocated(real_gdf_param))               deallocate(real_gdf_param)
if(allocated(real_weight_param))            deallocate(real_weight_param)
if(allocated(real_pdf_param))               deallocate(real_pdf_param)
if(allocated(jorek_filename))               deallocate(jorek_filename)
pdf_to_use => NULL(); weight_to_use => NULL();
gdf_to_use => NULL(); gdf_sampler_to_use => NULL();
call sim%finalize()
write(*,*) "Test: generalised_initialisation_gc_H_mu_psi completed."
!-------------------------------------------------------------------------------------------------

contains 
!> allocate a particle list using differeft particle types
!> inputs:
!>   particle_type: (integer) codify the particle type: 1-particle_gc (default)
!>                   2- particle_kinetic_leapfrog, 3- particle_gc_vpar
!>   n_particles:    (integer) number of particles
!>   particles:      (particle_base)(n_particles) particle list to be allocated
!> outputs:
!>   particles: (particle_base)(n_particles) allocated particle list
subroutine allocate_particle_list(particle_type,n_particles,particles)
use mod_particle_types, only: particle_base,particle_gc
use mod_particle_types, only: particle_kinetic_leapfrog,particle_gc_vpar
!> inputs:
integer,intent(in) :: particle_type,n_particles
!> inputs-outputs:
class(particle_base),dimension(:),allocatable :: particles
if(allocated(particles)) deallocate(particles)
if(particle_type.eq.1) then
  allocate(particle_gc::particles(n_particles))
elseif(particle_type.eq.2) then
  allocate(particle_kinetic_leapfrog::particles(n_particles))
elseif(particle_type.eq.3) then
  allocate(particle_gc_vpar::particles(n_particles))
else
  allocate(particle_gc::particles(n_particles))
endif
end subroutine allocate_particle_list

!> pdf to mimic, in this case just a dummy function returning 1
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   fields:       (fields_base) jorek MHD fields
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   pdf: (real8) value of the probability density, equal to 1
function pdf_psi_H_mu(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(pdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  pdf = 1d0
end function pdf_psi_H_mu

!> upper extremum of the pdf to mimic
!> inputs:
!>   nx:           (integer) number of variables
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   sup_pdf:      (real8) value of the probability density upper bound
function sup_pdf_psi_H_mu(nx,x_min,x_max,n_real_param,real_param,&
n_int_param,int_param) result(sup_pdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: sup_pdf
  sup_pdf = 1d0
end function sup_pdf_psi_H_mu

!> probability density function of the particle sampler,
!> in this case just a dummy function returning 1
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   fields:       (fields_base) jorek MHD fields
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution weight (normally mass/volume)
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   gdf: (real8) value of the sampler probability density 
function gdf_psi_H_mu(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(gdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: gdf
  gdf = 1d0
end function gdf_psi_H_mu

!> upper extremum of the sampler of the particle sampler
!> inputs:
!>   nx:           (integer) number of variables
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution weight (normally mass/volume)
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   sup_gdf:      (real8) value of the probability density upper bound
function sup_gdf_psi_H_mu(nx,x_min,x_max,n_real_param,real_param,&
n_int_param,int_param) result(sup_gdf)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: sup_gdf
  sup_gdf = 1d0
end function sup_gdf_psi_H_mu

!> The test code is used for emulating the initialise_particles_H_mu_psi
!> method using the initialise_particles_in_phase_space method
!> for guiding centers
!> inputs:
!>   n_x:          (integer) number of variables: must be 7
!>   x:            (real8)(n_x) uniform random numbers in [0,1)
!>   st:           (real8)(2) local particle coordinates 
!>   i_elm:        (integer) jorek element
!>   fields:       (fields_base) JOREK MHD field class
!>   x_min:        (real8)(n_x) lower bound uniform sampling, 1) poloidal flux
!>                 2) poloidal angle, 3) toroidal angle, 4) not used, 5) not used
!>                 6) gyro angle, 7) charge
!>   x_max:        (real8)(n_x) upper bound uniform sampling, 1) poloidal flux
!>                 2) poloidal angle, 3) toroidal angle, 4) not used, 5) not used
!>                 6) gyro angle, 7) charge
!>   n_real_param: (integer) size of the real parameters, must be: 2*n_elements+3
!>   real_param:   (real8)(2*n_elements+3) real gdf_sampler parameters:
!>                 1:n_elements              <- psi_minmax(:,1)
!>                 n_elements+1:2*n_elements <- psi_minmax(:,2)
!>                 2*n_elements+1            <- magnetic axis major radius
!>                 2*n_elements+2            <- magnetic axis vertical position
!>                 2*n_elements+3            <- particle mass
!>   n_int_param:  (integer) size of the integer parameters, must be: 1
!>   int_param:    (integer)(1) integer parameters: 1 <- if 1 paralle velocity is included
!> outputs:
!>   i_elm:        (integer) jorek element
!>   st:           (real8)(2) local particle coordinates 
!>   x:            (real8)(n_x) guiding center variables: 1) R, 2) Z, 3) phi,
!>                 4) energy, 5) magnetic moment, 6) gyro angle, 7) charge
subroutine gdf_psi_H_mu_sampler(n_x,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param)
  use constants,          only: MU_ZERO,EL_CHG,ATOMIC_MASS_UNIT
  use phys_module,        only: central_density
  use mod_model_settings, only: var_T,var_Vpar
  use mod_sampling,       only: sample_chi_squared_3
  use mod_fields,         only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in) :: n_x,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)            :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  class(fields_base),intent(in)               :: fields
  !> Inputs-Outputs:
  integer,intent(inout)               :: i_elm
  real*8,dimension(2),intent(inout)   :: st
  real*8,dimension(n_x),intent(inout) :: x
  !> Variables:
  real*8 :: psi,normB,electric_potential,u,temperature_ev
  real*8 :: R,R_s,R_t,Z,Z_s,Z_t
  real*8,dimension(1) :: P,P_s,P_t,P_phi
  real*8,dimension(3) :: B,E

  !> compute the psi,theta,phi coordinates
  x(1:3) = x_min(1:3)+(x_max(1:3)-x_min(1:3))*x(1:3)
  call find_theta_psi(fields%node_list,fields%element_list,[real_param(1:fields%element_list%n_elements),\
  real_param(fields%element_list%n_elements+1:2*fields%element_list%n_elements)],\
  x(2),x(1),x(3),real_param(2*fields%element_list%n_elements+1),\
  real_param(2*fields%element_list%n_elements+2),i_elm,st(1),st(2),R,Z); x(1:2) = [R,Z];
  if(i_elm.gt.0) then
    call fields%calc_EBpsiU(time,i_elm,st,x(3),E,B,psi,electric_potential); normB = norm2(B);
    !> compute the temperature
    call interp_PRZ(fields%node_list,fields%element_list,i_elm,[var_T],1,st(1),st(2),x(3),P,P_s,P_t,P_phi,&
    R,R_s,R_t,Z,Z_s,Z_t)
    temperature_ev = P(1)/(2d0*MU_ZERO*central_density*1d20*EL_CHG) !< temperature [eV]
#ifdef WITH_TiTe
    temperature_ev = 2d0*temperature_ev
#endif
    x(4) = temperature_ev*0.5d0*sample_chi_squared_3(x(4))
    u = 2*mod(x(5),0.5d0)
    x(5) = sign(x(4)*(2d0*u-u**2),x(5)-0.5d0)
    !> include parallel velocity if required
    if(int_param(1).eq.1) then
      call interp_PRZ(fields%node_list,fields%element_list,i_elm,[var_Vpar],1,st(1),st(2),x(3),P,P_s,P_t,P_phi,&
      R,R_s,R_t,Z,Z_s,Z_t)
      P(1) = P(1)*sqrt(real_param(2*fields%element_list%n_elements+3)*ATOMIC_MASS_UNIT/EL_CHG)
      x(4) = x(4) + P(1)*(P(1) + 2d0*sign(sqrt(x(4)-x(5)),x(5)))
    endif
      x(5) = x(5)/normB
      !> uniform distribution for both the gyro angle and the charge distribution
      !> very easy to modify for using the coronal charge equilibrium instead
      x(6:7) = x_min(6:7) + (x_max(6:7)-x_min(6:7))*x(6:7) 
  endif
end subroutine gdf_psi_H_mu_sampler

!> Compute the particle weight as in the H_mu_psi particle initialisation
!> basically, return 1
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   time:         (real8) physical time at wich the particle is sampled
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   fields:       (fields_base) jorek MHD fields
!>   n_real_param: (integer) N# of real input parameters
!>   real_param:   (real8)(n_real_param) real weight parameters
!>   n_int_param:  (integer) N# of integer input parameters
!>   int_param:    (integer)(n_int_param) integer weight parameters
!> outputs:
!>   weight:       (real) particle weight
function particle_weight_psi_H_mu(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param) result(weight)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: nx,i_elm,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(2),intent(in)              :: st
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: weight
  weight = 1d0
end function particle_weight_psi_H_mu

!> method used for copy a sample of the H mu psi sampler
!> into a partichel_gc particle type
!> the order of the variables in a sample are:
!> 1: R ,2: Z, 3: phi, 4: momentum, 5: pitch angle, 
!> 6: gyro-angle, 7: charge.
!> inputs: 
!>   p_inout:      (particle_base) particle to be initialised
!>   n_x:          (integer) size of the phase space sample: 7
!>   x:            (real8)(n_x) H mu psi sample
!>   time:         (real8) time of the simulation
!>   fields:       (fields_base) JOREK MHD fields
!>   n_real_param: (integer) number of real parameters: 1
!>   real_param:   (real8)(n_real_param) real parameters: 1: mass
!>   n_int_param:  (integer) number of integer parameters: 0
!>   int_param:    (integer)(n_real_param) integer parameters: empty
!> outputs:
!>   p_inout: (particle_base) initialised particle
subroutine copy_H_mu_psi_sample_to_p_gc(p_inout,&
n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
  use mod_particle_types,        only: particle_base
  use mod_particle_types,        only: particle_gc,particle_gc_vpar
  use mod_particle_types,        only: particle_kinetic_leapfrog
  use mod_boris,                 only: gc_to_kinetic_leapfrog
  use mod_gc_variational,        only: convert_gc_to_gc_vpar
  use mod_fields,                only: fields_base
  implicit none
  !> Inputs-Outputs:
  class(particle_base),intent(inout) :: p_inout
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: n_x,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)            :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> variables 
  type(particle_gc)   :: p_gc
  real*8              :: psi,electric_potential
  real*8,dimension(3) :: Evec,Bvec
  p_gc%x = p_inout%x; p_gc%st = p_inout%st; p_gc%i_elm = p_inout%i_elm;
  p_gc%weight = p_inout%weight; p_gc%i_life = p_inout%i_life;
  p_gc%t_birth = p_inout%t_birth;
  p_gc%E = x(4); p_gc%mu = x(5); p_gc%q = int(x(7),kind=1);
  select type (p=>p_inout)
  type is (particle_gc)
     p = p_gc
  type is (particle_kinetic_leapfrog)
    call fields%calc_EBpsiU(time,p%i_elm,p%st,p%x(3),Evec,Bvec,psi,electric_potential);
     p = gc_to_kinetic_leapfrog(p_gc,fields%node_list,fields%element_list,x(6),[0d0,0d0,0d0],Bvec,real_param(1),dt=0d0)
     p%q = p_gc%q
   type is (particle_gc_vpar)
     call fields%calc_EBpsiU(time,p%i_elm,p%st,p%x(3),Evec,Bvec,psi,electric_potential);
     call convert_gc_to_gc_vpar(p_gc,norm2(Bvec),real_param(1),p) 
     p%q = p_gc%q
  endselect
end subroutine copy_H_mu_psi_sample_to_p_gc

!> extract minimum and maximum poloidal flux for each mesh element
!> inputs:
!>   fields: (fields_base) JOREK MHD fields
!> outputs:
!>   psi_minmax_list_1d: (real8)(2*n_elements) mimimum (1:n_elements)
!>                       and maximum (n_elements+1:2*n_elements)
!>                       poloidal flux values for each element
subroutine extract_element_psi_minmax(fields,psi_minmax_list_1d)
  use mod_fields, only: fields_base
  implicit none
  !> inputs:
  class(fields_base),intent(in)  :: fields
  !> outputs:
  real*8,dimension(2*fields%element_list%n_elements),intent(out) :: psi_minmax_list_1d
  !> variables:
  integer :: ii
  real*8,dimension(fields%element_list%n_elements,2) :: psi_minmax_list
  !> initialisation
  psi_minmax_list(:,1) = 1d10; psi_minmax_list(:,2) = -1d10;
  !> extract maximum and minimun
  !$omp parallel do default(shared) private(ii)
  do ii=1,fields%element_list%n_elements
    call psi_minmax(fields%node_list,fields%element_list,ii,psi_minmax_list(ii,1),&
    psi_minmax_list(ii,2))
  enddo
  !$omp end parallel do
  psi_minmax_list_1d(1:fields%element_list%n_elements) = psi_minmax_list(:,1)
  psi_minmax_list_1d(fields%element_list%n_elements+1:2*fields%element_list%n_elements) = &
  psi_minmax_list(:,2)
end subroutine extract_element_psi_minmax

end program test_generalised_initialisation_gc_H_mu_psi
