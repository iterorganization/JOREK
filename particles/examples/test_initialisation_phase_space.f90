program test_initialisation_phase_space
!> program used for testing the initialisation method: 
!> initialise_particles_in_phase_space 
!> contained in mod_particle_group_init
!> using the relativistic kinetic particle model.
!> The idea is to provide a reject_uniform function
!> which reject a uniformly distributed random sample
!> between [0 1] as a function of a user defined 
!> 6D distribution. The distribution function used hereafter
!> is expressed in R,Z,phi coordinates for the physical space, 
!> momentum, pitch and gyro angles coordinates for the
!> momentum space and the cartesian coordinates for the 
!> distribution of electric charges. Be aware that 
!> the electric charge should be a double 
use constants,       only: TWOPI,PI
use phys_module,     only: xcase,xpoint
use data_structure,  only: type_bnd_node_list,type_bnd_element_list
use mod_boundary,    only: boundary_from_grid
use mod_expression,  only: exprs_all_int,init_expr,exprs,SI_UNITS
use mod_integrals3D, only: int3d_new
use mod_boundary,    only: boundary_from_grid
use equil_info
use mod_random_seed
use particle_tracer
use mod_particle_io, only: write_simulation_hdf5
implicit none

!> Variable declarations --------------------------------------------------------------------
type(pcg32_rng)              :: rng_pcg32
type(event)                  :: field_reader
type(type_bnd_node_list)     :: bnd_node_list
type(type_bnd_element_list)  :: bnd_elm_list
logical                      :: write_txt 
integer                      :: ii,n_variables,n_particles,nR,nZ,nphi,np,npitch,nchi
integer                      :: n_int_pdf_param,n_real_pdf_param,ifail
integer                      :: n_int_weight_param,n_real_weight_param
integer                      :: n_int_gdf_param,n_real_gdf_param
integer                      :: n_int_pdf_to_part_coord_param,n_real_pdf_to_part_coord_param
integer,dimension(:),allocatable :: int_pdf_param,int_weight_param,int_gdf_param
integer,dimension(:),allocatable :: int_pdf_to_part_coord_param
integer,dimension(:,:,:,:,:,:),allocatable :: histo
real*8                       :: start_time,mass,charge,error,error_norm
real*8                       :: error_avg_norm,pdf_upper_bound,gdf_upper_bound
real*8                       :: sup_pdf_safety_factor,n_tot_phys_particles
real*8                       :: error_n_phys_particles,error_n_phys_particles_norm
real*8,dimension(2)          :: Rbox,Zbox,Rbound,Zbound,Phibound
real*8,dimension(2)          :: Ekinbound,Pbound,Pitchbound,Chibound,Chargebound
real*8,dimension(7,2)        :: phase_space_bounds
real*8,dimension(:),allocatable           :: real_pdf_to_part_coord_param
real*8,dimension(:),allocatable           :: Rmesh,Zmesh,phimesh,pmesh,pitchmesh,chimesh
real*8,dimension(:),allocatable           :: real_pdf_param,real_weight_param,real_gdf_param
real*8,dimension(:),allocatable           :: DUMMY_REAL_ARRAY
real*8,dimension(:,:,:,:,:,:),allocatable :: expected_pdf,pdf_at_midpoints
character(len=125)                        :: test_case,particle_filename,mesh_filename_root
character(len=125)                        :: particle_pdf_filename,exact_pdf_filename
character(len=125)                        :: particle_histo_filename
character(len=15)                         :: particle_restart_filename
character(len=:),allocatable              :: jorek_filename
procedure(real_f),pointer                 :: pdf_to_use         => NULL()
procedure(real_f),pointer                 :: weight_to_use      => NULL()
procedure(real_f),pointer                 :: gdf_to_use         => NULL()
procedure(real_arr_inout_s),pointer       :: gdf_sampler_to_use => NULL()
procedure(part_inout_s),pointer           :: pdf_to_part_coord  => NULL()

!> MPI and groups initialisation ------------------------------------------------------------
call sim%initialize(num_groups=1)

!>-------------------------------------------------------------------------------------------
!> Define inputs ----------------------------------------------------------------------------
write_txt   = .false.
test_case   = 'jorek_current_density_re'
n_variables = 7
n_particles = 10000000
nR          = 2
nZ          = 2
nphi        = 2
np          = 81
npitch      = 81
nchi        = 81
start_time  = 0.d0
mass        = 5.48579909065d-4 !< electron mass in AMU
Rbound      = [0.d0,9.99d2]
Zbound      = [-9.99d2,9.99d2]
Phibound    = 2.5d-1*[TWOPI,5d0*PI]
Ekinbound   = [2d7-1d4,2d7+1d4]
Pitchbound  = [PI-2.95d-1,PI]
Chibound    = [0.d0,TWOPI]
Chargebound = -1.d0
charge      = -1.d0
allocate(character(len=25)::jorek_filename)
jorek_filename          = 'jorek_equilibrium' 
particle_filename       = 'jorek_particle_outputs.txt'
mesh_filename_root      = 'jorek_kinetic_mesh_'
particle_histo_filename = 'jorek_histogram_from_particles.txt'
particle_pdf_filename   = 'jorek_pdf_from_particles.txt'
exact_pdf_filename      = 'jorek_exact_pdf_at_midpoints.txt'
particle_restart_filename = 'part_restart.h5'
n_int_pdf_param         = 0
n_real_pdf_param        = 0
sup_pdf_safety_factor   = 1d0

!> Initialisation ---------------------------------------------------------------------------
!> allocate arrays and initialise them to 0
write(*,*) "Test: initialise_particle_in_phase_space: started."
write(*,*) "... setting-up test features"
allocate(Rmesh(nR)); allocate(Zmesh(nZ)); allocate(phimesh(nphi)); allocate(pmesh(np));
allocate(pitchmesh(npitch)); allocate(chimesh(nchi));
allocate(histo(nR,nZ,nphi,np,npitch,nchi));
allocate(expected_pdf(nR,nZ,nphi,np,npitch,nchi));
allocate(pdf_at_midpoints(nR,nZ,nphi,np,npitch,nchi));

Rmesh = 0.d0; Zmesh = 0.d0; phimesh = 0.d0; pmesh = 0.d0; pitchmesh = 0.d0; chimesh = 0.d0;
expected_pdf = 0.d0; pdf_at_midpoints = 0.d0;
!> read jorek field
field_reader = event(read_jorek_fields_interp_linear(basename=trim(jorek_filename),i=-1))
call with(sim,field_reader)
!> update equilibrium state
if(sim%my_id.eq.0) call boundary_from_grid(sim%fields%node_list,sim%fields%element_list,bnd_node_list,bnd_elm_list,.false.)
call broadcast_boundary(sim%my_id,bnd_elm_list,bnd_node_list)
call update_equil_state(sim%my_id,sim%fields%node_list,sim%fields%element_list,bnd_elm_list,xpoint,xcase)
!> initailise simulation parameters
sim%time = start_time
sim%groups(1)%mass = mass
allocate(particle_kinetic_relativistic::sim%groups(1)%particles(n_particles))
call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,Rbox(1),Rbox(2),Zbox(1),Zbox(2))
if(Rbox(1).ge.Rbound(1)) Rbound(1) = Rbox(1)
if((Rbox(2).lt.Rbound(2)).and.((Rbox(2)-Rbound(1)).gt.0.d0)) Rbound(2) = Rbox(2)
if((Zbox(1).gt.0.d0).and.(Zbound(1).ge.Zbox(1))) Zbound(1) = Zbox(1)
if((Zbox(1).lt.0.d0).and.(Zbound(1).lt.Zbox(1))) Zbound(1) = Zbox(1)
if((Zbox(2).gt.0.d0).and.(Zbound(2).ge.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
if((Zbox(2).lt.0.d0).and.(Zbound(2).lt.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
Pbound = mass*SPEED_OF_LIGHT*sqrt(((EL_CHG*Ekinbound/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT**2))+1.d0)**2-1.d0)
phase_space_bounds(:,1) = [Rbound(1),Zbound(1),Phibound(1),Pbound(1),Pitchbound(1),Chibound(1),charge]
phase_space_bounds(:,2) = [Rbound(2),Zbound(2),Phibound(2),Pbound(2),Pitchbound(2),Chibound(2),charge]
!> select the pdf to use
write(*,*) ' '
write(*,*) 'SELECT PDF TO USE: '
if(trim(test_case)=='jorek_current_density_re') then
  write(*,*) 'SELECTED: JOREK CURRENT DENSITY RUNAWAY ELECTRONS!'
  n_real_pdf_param = 3; allocate(real_pdf_param(n_real_pdf_param));
  real_pdf_param   = [1.d0,mass,sup_pdf_safety_factor]
  n_int_pdf_param  = 1; allocate(int_pdf_param(n_int_pdf_param));
  n_real_weight_param = 3; allocate(real_weight_param(n_real_weight_param));
  int_pdf_param(1)    = sim%my_id
  n_int_pdf_to_part_coord_param  = 0
  n_real_pdf_to_part_coord_param = 0
  pdf_to_use          => pdf_current_density_uniform_phase
  weight_to_use       => particle_weight_current_density_uniform_phase
  pdf_to_part_coord   => spherical_p_cartesian_q_to_relativistic_kinetic
  pdf_upper_bound     = sup_pdf_current_density_uniform_phase(n_variables,&
  phase_space_bounds(:,1),phase_space_bounds(:,2),sim%fields,&
  n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
  !> compute the integral of the current density in the volume
  allocate(DUMMY_REAL_ARRAY(n_real_weight_param+1)); call init_expr;
  call int3d_new(sim%my_id,sim%fields%node_list,sim%fields%element_list,bnd_node_list,bnd_elm_list,&
  exprs('int3d_jR_tot',1,exprs_all_int%n_coord,exprs_all_int),DUMMY_REAL_ARRAY,SI_UNITS)
  real_weight_param = [DUMMY_REAL_ARRAY(2),real(n_particles,kind=8),sim%groups(1)%mass]; 
  deallocate(DUMMY_REAL_ARRAY);
  call MPI_Bcast(real_weight_param,n_real_weight_param,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  !> compute the total number of physical particles
  n_tot_phys_particles = n_physical_particle_current_density_uniform_phase(n_variables,&
  start_time,sim%fields,phase_space_bounds(:,1),phase_space_bounds(:,2),&
  n_real_weight_param,real_weight_param,n_int_weight_param,int_weight_param)
elseif(trim(test_case)=='uniform_weight') then
  write(*,*) 'SELECTED: WEIGHTED PDF UNIFORM!'
  n_real_pdf_param = 2; allocate(real_pdf_param(n_real_pdf_param));
  n_int_pdf_to_part_coord_param  = 0
  n_real_pdf_to_part_coord_param = 0
  pdf_to_use        => pdf_uniform
  weight_to_use     => weight_uniform
  pdf_to_part_coord => spherical_p_cartesian_q_to_relativistic_kinetic
  pdf_upper_bound   = sup_pdf_uniform(n_variables-1,&
  phase_space_bounds(1:n_variables-1,1),phase_space_bounds(1:n_variables-1,2),&
  n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
  !> compute the plasma volume
  allocate(DUMMY_REAL_ARRAY(n_real_weight_param+1)); call init_expr;
  call int3d_new(sim%my_id,sim%fields%node_list,sim%fields%element_list,bnd_node_list,bnd_elm_list,&
  exprs('volume',1,exprs_all_int%n_coord,exprs_all_int),DUMMY_REAL_ARRAY,SI_UNITS)
  real_weight_param = [DUMMY_REAL_ARRAY(2),real(n_particles,kind=8)]; 
  deallocate(DUMMY_REAL_ARRAY);
  call MPI_Bcast(real_weight_param,n_real_weight_param,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  !> compute total number of physical particles
  n_tot_phys_particles = n_physical_particle_weight_uniform(n_variables,start_time,&
  sim%fields,phase_space_bounds(:,1),phase_space_bounds(:,2),&
  n_real_weight_param,real_weight_param,n_int_weight_param,int_weight_param)
else
  write(*,*) 'SELECTED: PDF UNIFORM (DEFAULT)!'
  n_int_pdf_to_part_coord_param  = 0
  n_real_pdf_to_part_coord_param = 0
  pdf_to_use        => pdf_uniform
  weight_to_use     => weight_uniform_one
  pdf_to_part_coord => spherical_p_cartesian_q_to_relativistic_kinetic
  pdf_upper_bound   = sup_pdf_uniform(6,phase_space_bounds(1:n_variables-1,1),&
  phase_space_bounds(1:n_variables-1,2),n_real_pdf_param,&
  real_pdf_param,n_int_pdf_param,int_pdf_param)
  !> compute total number of physical particles
  n_tot_phys_particles = real(n_particles,kind=8)
endif
!> select gdf to use
n_real_gdf_param = 0; n_int_gdf_param = 0;
gdf_to_use         => gdf_uniform_phase
gdf_sampler_to_use => gdf_uniform_sampler
gdf_upper_bound    = sup_gdf_uniform_phase(n_variables-1, &
phase_space_bounds(1:n_variables-1,1),phase_space_bounds(1:n_variables-1,2), &
n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
write(*,*) ' '

!> Test particle initialisation -------------------------------------------------------------
write(*,*) "... initialising particles in phase space"
call initialise_particles_in_phase_space(n_variables,sim%groups(1)%particles,sim%fields,rng_pcg32,&
pdf_to_use,weight_to_use,gdf_to_use,gdf_sampler_to_use,pdf_upper_bound,gdf_upper_bound,&
pdf_to_part_coord,sim%groups(1)%mass,start_time,phase_space_bounds,n_real_pdf_param,real_pdf_param,&
n_int_pdf_param,int_pdf_param,n_real_weight_param,real_weight_param,n_int_weight_param,&
int_weight_param,n_real_gdf_param,real_gdf_param,n_int_gdf_param,int_gdf_param,&
n_real_pdf_to_part_coord_param,real_pdf_to_part_coord_param,&
n_int_pdf_to_part_coord_param,int_pdf_to_part_coord_param)

!> Produce the expected pdf fromt the particle histogram --------------------------------------
write(*,*) "... building particle histogram and computing the expected pdf"
select type(plist=>sim%groups(1)%particles)
  type is (particle_kinetic_relativistic)
  call compute_cylindrical_spherical_histogram_pdf(histo,expected_pdf,Rmesh,&
  Zmesh,phimesh,pmesh,pitchmesh,chimesh,n_particles,plist,start_time,nR,nZ,nphi,np,&
  npitch,nchi,Rbound,Zbound,Phibound,Pbound,Pitchbound,Chibound,sim%fields)
end select

!> Compute the input pdf at the mesh element midpoints and the L2 error w.r.t. the 
!> the expected pdf from the  particle histogram and the error on the 
!> total number of physical particles
write(*,*) "... computing the input pdf at the midpoints of the mesh elements"
call evaluate_pdf_at_midpoints(pdf_at_midpoints,nR,nZ,nphi,np,npitch,nchi,Rmesh,Zmesh,phimesh,&
     pmesh,pitchmesh,chimesh,charge,start_time,pdf_to_use,sim%fields,n_real_pdf_param,&
     real_pdf_param,n_int_pdf_param,int_pdf_param)
write(*,*) "... computing L2 error"
call compute_error_norm2_ndim6(error,error_norm,error_avg_norm,nR-1,nZ-1,nphi-1,np-1,&
npitch-1,nchi-1,expected_pdf,pdf_at_midpoints,pdf_upper_bound) 
call compute_error_tot_n_phys_particles(error_n_phys_particles,error_n_phys_particles_norm,&
n_tot_phys_particles,sim%groups(1)%particles)

!> Log test results -------------------------------------------------------------------------
write(*,*) "... logging test results"
write(*,*) "L2 error between the expected pdf from particle histogram and the input pdf at mid points: ",error
write(*,*) "L2 error normalized to the maximum of the input pdf: ",error_norm
write(*,*) "L2 error averaged and normalised to the maximum of the input pdf: ",error_avg_norm
write(*,*) "Error on the total number of physical particles: ",error_n_phys_particles
write(*,*) "Error on the total number of physical particles normalised: ",error_n_phys_particles_norm
write(*,*) " "

!> Write data into file ---------------------------------------------------------------------
write(*,*) "... write particle restart file"
call write_simulation_hdf5(sim,trim(particle_restart_filename))
write(*,*) "... writing data in files"
if(write_txt) then
  call dump_kinetic_particles_in_txt(n_particles,sim%groups(1)%particles,sim%groups(1)%mass,&
  start_time,sim%fields,particle_filename,ifail)
  call write_array1d_double(nR,Rmesh,trim(trim(mesh_filename_root)//'R.txt'),ifail)
  call write_array1d_double(nZ,Zmesh,trim(trim(mesh_filename_root)//'Z.txt'),ifail)
  call write_array1d_double(nphi,phimesh,trim(trim(mesh_filename_root)//'phi.txt'),ifail)
  call write_array1d_double(np,pmesh,trim(trim(mesh_filename_root)//'p.txt'),ifail)
  call write_array1d_double(npitch,pitchmesh,trim(trim(mesh_filename_root)//'pitch.txt'),ifail)
  call write_array1d_double(nchi,chimesh,trim(trim(mesh_filename_root)//'gyro.txt'),ifail)
  call write_array6d_integer(nchi-1,npitch-1,np-1,nphi-1,nZ-1,nR-1,histo,trim(particle_histo_filename),ifail)
  call write_array6d_double(nchi-1,npitch-1,np-1,nphi-1,nZ-1,nR-1,expected_pdf,trim(particle_pdf_filename),ifail)
  call write_array6d_double(nchi-1,npitch-1,np-1,nphi-1,nZ-1,nR-1,pdf_at_midpoints,trim(exact_pdf_filename),ifail)
endif
!> add here method for writing the charge mesh distribution as well
!> Clean-up ---------------------------------------------------------------------------------
deallocate(Rmesh); deallocate(Zmesh); deallocate(phimesh); deallocate(pmesh); 
deallocate(pitchmesh); deallocate(chimesh); deallocate(histo); deallocate(expected_pdf); 
deallocate(pdf_at_midpoints); if(allocated(real_pdf_param)) deallocate(real_pdf_param);
if(allocated(int_pdf_param)) deallocate(int_pdf_param);
if(allocated(int_weight_param)) deallocate(int_weight_param);
if(allocated(real_weight_param)) deallocate(real_weight_param);
pdf_to_use => NULL(); weight_to_use => NULL();
gdf_to_use => NULL(); gdf_sampler_to_use => NULL();
call sim%finalize()
write(*,*) "Test: initialise_particle_in_phase_space: completed."

contains

!> Compute the L2 error of 6D-arrays
subroutine compute_error_norm2_ndim6(error_L2,error_L2_norm,error_L2_avg_norm,&
n1,n2,n3,n4,n5,n6,array1,array2,sup_array2)
  implicit none
  !> inputs
  integer,intent(in)                             :: n1,n2,n3,n4,n5,n6
  real*8,intent(in)                              :: sup_array2
  real*8,dimension(n1,n2,n3,n4,n5,n6),intent(in) :: array1,array2
  !> outputs
  real*8,intent(out) :: error_L2,error_L2_norm,error_L2_avg_norm
  !> variables
  integer :: ii,jj,kk,pp,qq
  real*8  :: error

  !> initialisation
  error = 0.d0
  !> compute norm2 error, the openmp collapse clause is used at first. Manual collapse
  !> of the indices should be tried in case of reduced performance.
  !$omp parallel do default(none) firstprivate(n2,n3,n4,n5,n6) &
  !$omp shared(array1,array2) private(ii,jj,kk,pp,qq) reduction(+:error) &
  !$omp collapse(5)
  do ii=1,n6
    do jj=1,n5
      do kk=1,n4
        do pp=1,n3
          do qq=1,n2
            error = error + dot_product((array2(:,qq,pp,kk,jj,ii)-array1(:,qq,pp,kk,jj,ii)),&
            (array2(:,qq,pp,kk,jj,ii)-array1(:,qq,pp,kk,jj,ii)))
          enddo
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do
  error_L2 = sqrt(error); error_L2_norm = error_L2/abs(sup_array2);
  error_L2_norm = sqrt(error/(n1*n2*n3*n4*n5*n6))/abs(sup_array2)
end subroutine compute_error_norm2_ndim6

!> compute the error between the expected total number of physical particles
!> from pdf integration and the actual total number of physical particles 
!> stored in the particles data structure
!> inputs:
!>   n_tot_phys_particles: (real8) expected number of physical particles
!>   particles:            (particle_base) the list of particles 
!>                         to be tested
!> outputs:
!>   error_n_phys_particles: (real8) error of the total number of physical particles
subroutine compute_error_tot_n_phys_particles(error_n_phys_particles,&
error_n_phys_particles_norm,n_tot_phys_particles,particles)
  use mod_particle_types, only: particle_base
  implicit none
  !> Inputs:
  real*8,intent(in)  :: n_tot_phys_particles
  class(particle_base),dimension(:),allocatable :: particles
  !> Outputs:
  real*8,intent(out) :: error_n_phys_particles,error_n_phys_particles_norm
  !> Variables:
  integer :: ii,n_num_particles
  real*8 :: n_phys_particles_num
  !> compute the total number of physical particles from 
  !> particle simulations
  n_num_particles = size(particles); n_phys_particles_num = 0.d0;
  !$omp parallel do default(none) shared(particles) private(ii) &
  !$omp firstprivate(n_num_particles) reduction(+:n_phys_particles_num)
  do ii=1,n_num_particles
    n_phys_particles_num = n_phys_particles_num + particles(ii)%weight
  enddo
  !$omp end parallel do
  error_n_phys_particles = abs(n_phys_particles_num-n_tot_phys_particles)
  error_n_phys_particles_norm = error_n_phys_particles/n_tot_phys_particles
end subroutine compute_error_tot_n_phys_particles

!> compute the particle histogram and estimate the pdf using bins of constant volume
!> inputs:
!>   n_particles: (integer) number of particles
!>   particles:   (particle_kinetic_relativistic)(n_particles) particle list
!>   time:        (real8) time of the required MHD fields
!>   nR:          (integer) size of the mesh along the major radius
!>   nZ:          (integer) size of the mesh along the vertical coordinate  
!>   nphi:        (integer) size of the mesh along the toroidal angle
!>   np:          (integer) size of the mesh along the total momentum
!>   npitch:      (integer) size of the mesh along the pitch angle
!>   ngyro:       (integer) size of the mesh along the gyro angle
!>   Rbound:      (real8)(2) upper and lower major radius bounds
!>   Zbound:      (real8)(2) upper and lower vertical coordinate bounds
!>   phibound:    (real8)(2) upper and lower toroidal angle bounds
!>   pbound:      (real8)(2) upper and lower total momentum bounds
!>   pitchbound:  (real8)(2) upper and lower pitch angle bounds
!>   gyrobound:   (real8)(2) upper and lower gyro angle bounds
!>   fields:      (fields_base) jorek MHD fields
!> outputs:
!>   histo:         (integer)(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1) 6D particle histogram
!>   estimated_pdf: (real8)(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1) 6D pdf estimation   
!>   Rmesh:         (real8)(nR) major radius mesh equidistant in R**2
!>   Zmesh:         (real8)(nZ) vertical coordinate equidistant mesh
!>   phimesh:       (real8)(nphi) toroidal angle equidistant mesh
!>   pmesh:         (real8)(np) total momentum mesh equidistant in p**3
!>   pitchmesh:     (real8)(npitch) pitch angle mesh equidistant in cos(pitch)
!>   gyromesh:      (real8)(ngyro) gyro angle equidistant mesh
subroutine compute_cylindrical_spherical_histogram_pdf(histo,estimated_pdf,Rmesh,&
Zmesh,phimesh,pmesh,pitchmesh,gyromesh,n_particles,particles,time,nR,nZ,nphi,np,&
npitch,ngyro,Rbound,Zbound,phibound,pbound,pitchbound,gyrobound,fields)
  use constants,                 only: TWOPI
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  use mod_fields,                only: fields_base
  use mod_particle_types,        only: particle_kinetic_relativistic
  use mod_pusher_tools,          only: get_orthonormals
  implicit none
  !> inputs:
  class(fields_base),intent(in) :: fields
  integer,intent(in) :: n_particles,nR,nZ,nphi,np,npitch,ngyro
  real*8,intent(in)  :: time
  real*8,dimension(2),intent(in) :: Rbound,Zbound,phibound
  real*8,dimension(2),intent(in) :: pbound,pitchbound,gyrobound
  type(particle_kinetic_relativistic),dimension(n_particles),intent(in) :: particles
  !> outputs:
  integer,dimension(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1),intent(out) :: histo
  real*8,dimension(nR),intent(out)     :: Rmesh
  real*8,dimension(nZ),intent(out)     :: Zmesh
  real*8,dimension(nphi),intent(out)   :: phimesh
  real*8,dimension(np),intent(out)     :: pmesh
  real*8,dimension(npitch),intent(out) :: pitchmesh
  real*8,dimension(ngyro),intent(out)  :: gyromesh 
  real*8,dimension(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1),intent(out) :: estimated_pdf
  !> variables
  integer              :: ii,jj
  integer,dimension(6) :: ids
  real*8               :: dvolume_real,dvolume_momentum,psi,U,ppar
  real*8,dimension(3)  :: pcyl,B,E,e1,e2
  real*8,dimension(6)  :: dist,coord,mesh_init
  !> Initialisation
  ids = 0; dist = 1d0; histo=0; estimated_pdf = 0;
  mesh_init = [Rbound(1)**2,Zbound(1),phibound(1),pbound(1)**3,&
  cos(pitchbound(1)),gyrobound(1)]
  !> Compute cylindrical mesh real space keeping the volume of each bin constant
  dvolume_real = (5d-1*(Rbound(2)**2-Rbound(1)**2)*&
  (Zbound(2)-Zbound(1))*(phibound(2)-phibound(1)))/&
  ((nR-1)*(nZ-1)*(nphi-1))
  call compute_equidistant_mesh(Zmesh,nZ,Zbound)
  call compute_equidistant_mesh(phimesh,nphi,phibound)
  dist(2) = Zmesh(2)-Zmesh(1); dist(3) = phimesh(2)-phimesh(1);
  dist(1) = ((2d0*dvolume_real)/(dist(2)*dist(3)));
  Rmesh(1) = Rbound(1)**2
  do ii=1,nR-1
    Rmesh(ii+1) = Rmesh(1) + real(ii,kind=8)*dist(1)
  enddo
  Rmesh = sqrt(Rmesh)
  !> Compute spherical mesh momentum space keeping the volume of each bin constant
  dvolume_momentum = ((pbound(2)**3-pbound(1)**3)*&
  (cos(pitchbound(1))-cos(pitchbound(2)))*(gyrobound(2)-gyrobound(1)))/&
  (3d0*(npitch-1)*(np-1)*(ngyro-1)) 
  call compute_equidistant_mesh(gyromesh,ngyro,gyrobound)
  call compute_equidistant_mesh(pitchmesh,npitch,cos(pitchbound))
  dist(5) = pitchmesh(2) - pitchmesh(1); dist(6) = gyromesh(2)-gyromesh(1);
  dist(4) = ((-3d0*dvolume_momentum)/(dist(5)*dist(6))) 
  pmesh(1) = pbound(1)**3
  do ii=1,np-1
    pmesh(ii+1) = pmesh(1) + real(ii,kind=8)*dist(4)
  enddo
  pmesh = pmesh**(1d0/3d0); pitchmesh = acos(pitchmesh);
  !> compute histogram
  !$omp parallel do default(shared) &
  !$omp firstprivate(n_particles,time,dist,mesh_init,npitch) &
  !$omp private(ii,ids,coord,E,B,psi,U,pcyl,e1,e2,ppar) &
  !$omp reduction(+:histo,estimated_pdf)
  do ii=1,n_particles
    ids = 1; coord(1:3) = particles(ii)%x
    !> compute the momentum in spherical coordinates
    call fields%calc_EBpsiU(time,particles(ii)%i_elm,particles(ii)%st,&
    coord(3),E,B,psi,U)
    pcyl = vector_cartesian_to_cylindrical(coord(3),particles(ii)%p)
    B = B/norm2(B)
    call get_orthonormals(B,e1,e2)
    coord(4)  = norm2(pcyl)
    ppar  = dot_product(pcyl,B)
    coord(5) = ppar/coord(4)
    coord(6)  = atan2(dot_product(pcyl-ppar*B,e2),dot_product(pcyl-ppar*B,e1))
    if(coord(6).lt.0.d0) coord(6) = TWOPI+coord(6)
    coord(6) = mod(coord(6),TWOPI)
    !> compute indices
    coord(1) = coord(1)**2; coord(4) = coord(4)**3;
    ids = ids+floor((coord-mesh_init)/dist)
    !> compute histogram and pdf
    histo(ids(1),ids(2),ids(3),ids(4),ids(5),ids(6)) = &
    histo(ids(1),ids(2),ids(3),ids(4),ids(5),ids(6)) + 1
    estimated_pdf(ids(1),ids(2),ids(3),ids(4),ids(5),ids(6)) = &
    estimated_pdf(ids(1),ids(2),ids(3),ids(4),ids(5),ids(6)) + particles(ii)%weight
  enddo 
  !$omp end parallel do
  estimated_pdf = estimated_pdf/(dvolume_real*dvolume_momentum)
end subroutine compute_cylindrical_spherical_histogram_pdf

!> Generate an equidistant mesh
!> inputs:
!>   n_points: (integer) number of nodes
!>   bounds:   (real8)(2) upper and lower mesh extrema
!> outputs:
!>   mesh:     (real8)(n_points) equidistant mesh nodes
subroutine compute_equidistant_mesh(mesh,n_points,bounds)
  implicit none
  !> input variables
  integer,intent(in)                     :: n_points
  real*8,dimension(2),intent(in)         :: bounds
  !> output variables
  real*8,dimension(n_points),intent(out) :: mesh
  !> variables
  real*8                                 :: delta
  
  delta = (bounds(2)-bounds(1))/(n_points-1)
  do ii=1,n_points
    mesh(ii) = bounds(1)+real(delta*(ii-1),kind=8)
  enddo
end subroutine compute_equidistant_mesh

!> evaluate the pdf at the mesh element midpoints
!> inputs:
!>   nR:        (integer) number of the major radius mesh nodes
!>   nZ:        (integer) number of the vertical mesh nodes
!>   nphi:      (integer) number of the toroidal angle mesh nodes
!>   np:        (integer) number of the momentum mesh nodes
!>   npitch:    (integer) number of the pitch-angle mesh nodes
!>   ngyro:     (integer) number of the gyro-angle mesh nodes
!>   Rmesh:     (real8)(nR) major radius mesh
!>   Zmesh:     (real8)(nZ) vertical coordinate mesh
!>   phimesh:   (real8)(nphi) toroidal mesh
!>   pitchmesh: (real8)(npitch) pitch angle mesh
!>   gyromesh:  (real8)(ngyro) gyro angle mesh
!>   pdf:       (procedure) particle distribution function
!>   fields:    (fields_base) jorek MHD fields
!> outputs:
!>   pdf_at_midpoints: (real8) value of the pdf at the midpoints
subroutine evaluate_pdf_at_midpoints(&
pdf_midpoints,nR,nZ,nphi,np,npitch,ngyro,Rmesh,Zmesh,phimesh,pmesh,& 
pitchmesh,gyromesh,charge,time,pdf,fields,n_real_pdf_param_in,&
real_pdf_param_in,n_int_pdf_param_in,int_pdf_param_in)
   use mod_fields, only: fields_base
  implicit none
  !> Parameters
  integer,parameter                   :: nx=7
  !> Inputs
  class(fields_base),intent(in)       :: fields
  integer,intent(in)                  :: nR,nZ,nphi,np,npitch,ngyro
  real*8,intent(in)                   :: charge,time
  real*8,dimension(nR),intent(in)     :: Rmesh
  real*8,dimension(nZ),intent(in)     :: Zmesh
  real*8,dimension(nphi),intent(in)   :: phimesh
  real*8,dimension(np),intent(in)     :: pmesh
  real*8,dimension(npitch),intent(in) :: pitchmesh
  real*8,dimension(ngyro),intent(in)  :: gyromesh
  integer,intent(in),optional         :: n_real_pdf_param_in,n_int_pdf_param_in
  integer,dimension(:),allocatable,intent(in),optional :: int_pdf_param_in
  real*8,dimension(:),allocatable,intent(in),optional  :: real_pdf_param_in
  procedure(real_f)                    :: pdf
  !> Outputs:
  real*8,dimension(nR-1,nZ-1,nphi-1,np-1,npitch-1,ngyro-1),intent(out) :: pdf_midpoints
  !> Variables:
  integer                  :: ii,jj,kk,pp,qq,rr
  integer                  :: i_elm,ifail,n_real_pdf_param,n_int_pdf_param
  integer,dimension(:),allocatable :: int_pdf_param
  real*8                   :: dummy_double_1,dummy_double_2
  real*8,dimension(2)      :: st
  real*8,dimension(nx)     :: x_midpoints,x_min,x_max
  real*8,dimension(:),allocatable :: real_pdf_param
  !> initialisation
  pdf_midpoints = 0.d0; x_midpoints = [0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,charge];
  x_min = [Rmesh(1),Zmesh(1),phimesh(1),pmesh(1),pitchmesh(1),gyromesh(1),charge]
  x_max = [Rmesh(nR),Zmesh(nZ),phimesh(nphi),pmesh(np),pitchmesh(npitch),gyromesh(ngyro),charge]
  n_real_pdf_param = 0; if(present(n_real_pdf_param_in)) n_real_pdf_param = n_real_pdf_param_in;
  if(present(real_pdf_param_in).and.(n_real_pdf_param_in.gt.0)) then
    allocate(real_pdf_param(n_real_pdf_param)); real_pdf_param = real_pdf_param_in;
  endif
  n_int_pdf_param = 0; if(present(n_int_pdf_param_in)) n_int_pdf_param = n_int_pdf_param_in;
  if(present(int_pdf_param_in).and.(n_int_pdf_param_in.gt.0)) then
    allocate(int_pdf_param(n_int_pdf_param)); int_pdf_param = int_pdf_param_in;
  endif
  !> loop over all midpoints, try with openmp collapse clause first. If slow,
  !> manually collapse all loops in one
  !$omp parallel do default(shared) &
  !$omp firstprivate(ngyro,npitch,np,nphi,nZ,nR,time,x_midpoints,x_min,x_max,&
  !$omp n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param) &
  !$omp private(ii,jj,kk,pp,qq,rr,i_elm,st,dummy_double_1,dummy_double_2) &
  !$omp collapse(6)
  do ii=1,ngyro-1
    do jj=1,npitch-1
      do kk=1,np-1
        do pp=1,nphi-1
          do qq=1,nZ-1
            do rr=1,nR-1
              !> compute midpoint in the local mesh cell
              x_midpoints(1:6) = 5.d-1*([Rmesh(rr+1),Zmesh(qq+1),&
              phimesh(pp+1),pmesh(kk+1),pitchmesh(jj+1),gyromesh(ii+1)] + &
              [Rmesh(rr),Zmesh(qq),phimesh(pp),pmesh(kk),pitchmesh(jj),gyromesh(ii)])
              !> find local JOREK mesh coordinates 
              call find_RZ(fields%node_list,fields%element_list,x_midpoints(1),&
              x_midpoints(2),dummy_double_1,dummy_double_2,i_elm,st(1),st(2),ifail)
              !> estimate the number of particles at the kinetic mesh midpoint
              pdf_midpoints(rr,qq,pp,kk,jj,ii) = pdf(nx,x_midpoints,st,time,i_elm,fields,&
              x_min,x_max,n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !$omp end parallel do
  !> Clean-up
  if(allocated(real_pdf_param)) deallocate(real_pdf_param);
  if(allocated(int_pdf_param))  deallocate(int_pdf_param);
end subroutine evaluate_pdf_at_midpoints

!> Dump particle kinetic list in txt file
subroutine dump_kinetic_particles_in_txt(n_particles,particles,mass,time,fields,filename,ifail)
  use constants,                 only: TWOPI
  use mod_particle_types,        only: particle_base
  use mod_particle_types,        only: particle_kinetic_relativistic
  use mod_fields,                only: fields_base
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  use mod_pusher_tools,          only: get_orthonormals
  implicit none
  !> Onput-Outputs
  integer,intent(inout) :: ifail
  !> Inputs
  integer,intent(in)                                     :: n_particles
  class(particle_base),dimension(n_particles),intent(in) :: particles
  class(fields_base),intent(in)                          :: fields
  character(len=*),intent(in)                            :: filename
  real*8,intent(in)                                      :: mass,time
  !> Variables
  integer :: ii
  real*8  :: ptot,ppar,pitch,gyro,psi,U
  real*8,dimension(3) :: B,E,pcyl,e1,e2

  !> open the file and write particle properties
  open(unit=42,file=trim(filename),action='write',blank='NULL',form='formatted',status='unknown',iostat=ifail)
  select type (plist=>particles)
  type is (particle_kinetic_relativistic)
    do ii=1,n_particles
    !> compute the momentum in spherical coordinates
    call fields%calc_EBpsiU(time,plist(ii)%i_elm,plist(ii)%st,plist(ii)%x(3),E,B,psi,U)
    pcyl = vector_cartesian_to_cylindrical(plist(ii)%x(3),plist(ii)%p)
    B = B/norm2(B)
    call get_orthonormals(B,e1,e2)
    ptot  = norm2(pcyl)
    ppar  = dot_product(pcyl,B)
    pitch = acos(ppar/ptot)
    gyro  = atan2(dot_product(pcyl-ppar*B,e2),dot_product(pcyl-ppar*B,e1))
    if(gyro.lt.0.d0) gyro = TWOPI+gyro
    gyro = mod(gyro,TWOPI)
      write(42,'(12E40.16E4)') plist(ii)%x(1),plist(ii)%x(2),plist(ii)%x(3),plist(ii)%p(1),&
      plist(ii)%p(2),plist(ii)%p(3),plist(ii)%weight,real(plist(ii)%q,kind=8),&
      mass,ptot,pitch,gyro
    enddo
  end select
  close(42) !< close file
end subroutine dump_kinetic_particles_in_txt

!> write 1d array of doubles in txt file
subroutine write_array1d_double(size1,array1d,filename,ifail)
  implicit none
  !> Inputs-Outputs
  integer,intent(inout) :: ifail
  !> Inputs:
  integer,intent(in) :: size1
  real*8,dimension(size1),intent(in) :: array1d
  character(len=*),intent(in) :: filename
  !> Variables
  integer :: ii
  !> write 1d array of doubles in txt file
  open(unit=43,file=trim(filename),action='write',blank='NULL',form='formatted',status='unknown',iostat=ifail) 
  do ii=1,size1
    write(43,'(E40.16e4)') array1d(ii)
  enddo
  close(43)
end subroutine write_array1d_double

!> write 6d array of doubles in txt file
subroutine write_array6d_double(size1,size2,size3,size4,size5,size6,array6d,filename,ifail)
  implicit none
  !> Inputs-Outputs
  integer,intent(inout) :: ifail
  !> Inputs
  integer,intent(in) :: size1,size2,size3,size4,size5,size6
  real*8,dimension(size1,size2,size3,size4,size5,size6),intent(in) :: array6d
  character(len=*),intent(in) :: filename
  !> Variables
  integer :: ii,jj,kk,pp,qq,rr
  !> write 6d array of doubles in array
  open(unit=44,file=trim(filename),action='write',blank='NULL',form='formatted',status='unknown',iostat=ifail)
  do rr=1,size6
    do qq=1,size5
      do pp=1,size4
        do kk=1,size3
          do jj=1,size2
            do ii=1,size1
              write(44,'(E40.16e4)') array6d(ii,jj,kk,pp,qq,rr)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  close(44)
end subroutine write_array6d_double

!> write 6d array of integer in txt file
subroutine write_array6d_integer(size1,size2,size3,size4,size5,size6,array6d,filename,ifail)
  implicit none
  !> Inputs-Outputs
  integer,intent(inout) :: ifail
  !> Inputs
  integer,intent(in) :: size1,size2,size3,size4,size5,size6
  integer,dimension(size1,size2,size3,size4,size5,size6),intent(in) :: array6d
  character(len=*),intent(in) :: filename
  !> Variables
  integer :: ii,jj,kk,pp,qq,rr
  !> write 6d array of doubles in array
  open(unit=45,file=trim(filename),action='write',blank='NULL',form='formatted',status='unknown',iostat=ifail)
  do rr=1,size6
    do qq=1,size5
      do pp=1,size4
        do kk=1,size3
          do jj=1,size2
            do ii=1,size1
              write(45,'(I6)') array6d(ii,jj,kk,pp,qq,rr)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  close(45)
end subroutine write_array6d_integer

!> Definitions of the probability density function (PDF) ------------------------ 

!> method used for transforming the momentum space from spherical
!> coordinates (p,pitch,chi) and cartesian coordinates for the
!> charge state into particle kinetic relativistic coordinates
!> the order of the variables in a sample are:
!> 1: R ,2: Z, 3: phi, 4: momentum, 5: pitch angle, 
!> 6: gyro angle, 7: charge.
!> inputs: 
!>   p_inout:      (particle_base) particle to be initialised
!>   n_x:          (integer) size of the phase space sample
!>   x:            (real8)(n_x) phase space sample in spehrical
!>                 momentum coordinates cartesian charge coordinates
!>   time:         (real8) time of the simulation
!>   fields:       (fields_base) JOREK MHD fields
!>   n_real_param: (integer) number of real parameters: 0
!>   real_param:   (real8)(n_real_param) real parameters: empty
!>   n_int_param:  (integer) number of integer parameters: 0
!>   int_param:    (integer)(n_real_param) integer parameters: empty
!> outputs:
!>   p_inout: (particle_base) initialised particle
subroutine spherical_p_cartesian_q_to_relativistic_kinetic(p_inout,&
n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
  use mod_particle_types,        only: particle_base
  use mod_particle_types,        only: particle_kinetic_relativistic
  use mod_fields,                only: fields_base
  use mod_pusher_tools,          only: get_orthonormals
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
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
  real*8              :: psi,U
  real*8,dimension(3) :: B_field,E_field,e1,e2
  select type (p=>p_inout)
  type is (particle_kinetic_relativistic)
    call fields%calc_EBpsiU(time,p%i_elm,p%st,p%x(3),&
    E_field,B_field,psi,U); B_field = B_field/norm2(B_field);
    call get_orthonormals(B_field,e1,e2)
    p%p = x(4)*(cos(x(5))*B_field + sin(x(5))*(cos(x(6))*e1 + sin(x(6))*e2))
    p%p = vector_cylindrical_to_cartesian(p%x(3),p%p)
    p%q = int(x(7),kind=1)
  end select
end subroutine spherical_p_cartesian_q_to_relativistic_kinetic

!> Phase space distribution for testing
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
!>   pdf: (real8) value of the probability density 
function pdf_uniform(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param)
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
  real*8 :: pdf_uniform
  !> Evalutate pdf
  pdf_uniform = 6.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))*&
  (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)*&
  (cos(x_min(5))-cos(x_max(5)))*(x_max(6)-x_min(6)))
  if(n_real_param.gt.0) pdf_uniform = real_param(1)*pdf_uniform
end function pdf_uniform

!> Upper bound of the phase space distribution for testing
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
!>   sup_pdf:      (real8) value of the probability density upper bound
function sup_pdf_uniform(nx,x_min,x_max,n_real_param,real_param,&
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
  !> Evalutate the upper extremum of the pdf
  sup_pdf = 6.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))*&
  (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)*&
  (cos(x_min(5))-cos(x_max(5)))*(x_max(6)-x_min(6)))
  if(n_real_param.gt.0) sup_pdf = real_param(1)*sup_pdf
end function sup_pdf_uniform

!> Phase space distribution based on the plasma current density
!> and uniform phase space distribution. The momentum distribution
!> is considered uniform for relativistic particle hence, the 
!> appearance of the relativistic factor.
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   time:         (real8) physical time at wich the particle is sampled
!>   fields:       (fields_base) jorek MHD fields
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution coefficient (not used)
!>                 2) particle mass in AMU
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   pdf: (real8) value of the probability density 
function pdf_current_density_uniform_phase(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param) result(pdf)
  use constants,          only: MU_ZERO,EL_CHG,SPEED_OF_LIGHT
  use mod_model_settings, only: var_zj
  use mod_interp,         only: interp_PRZ
  use mod_fields,         only: fields_base
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
  !> Variables:
  real*8 :: DUMMY_DOUBLE_1,DUMMY_DOUBLE_2
  real*8,dimension(1) :: jphi
  !> interpolate the jorek toroidal current density at the particle position
  call interp_PRZ(fields%node_list,fields%element_list,i_elm,[var_zj],1,&
  st(1),st(2),x(3),jphi,DUMMY_DOUBLE_1,DUMMY_DOUBLE_2)
  !> compute the pdf
  DUMMY_DOUBLE_1 = sqrt((x_max(4)**2)/((real_param(2)*SPEED_OF_LIGHT)**2)+1.d0)
  DUMMY_DOUBLE_2 = sqrt((x_min(4)**2)/((real_param(2)*SPEED_OF_LIGHT)**2)+1.d0)
  pdf = ((DUMMY_DOUBLE_1**3)-3.d0*DUMMY_DOUBLE_1) - ((DUMMY_DOUBLE_2**3)-3.d0*DUMMY_DOUBLE_2);
  pdf = pdf*(cos(x_min(5))**2 - cos(x_max(5))**2)*(x_max(6)-x_min(6))
  pdf =(-6.d0*jphi(1))/(pdf*x(7)*EL_CHG*MU_ZERO*(real_param(2)**3)*(SPEED_OF_LIGHT**4)*x(1))
end function pdf_current_density_uniform_phase

!> Upper bound phase space distribution based on the plasma current density
!> and uniform phase space distribution. The momentum distribution
!> is considered uniform for relativistic particle hence, the 
!> appearance of the relativistic factor.
!> inputs:
!>   nx:           (integer) number of variables
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   fields:       (fields_base) jorek MHD fields
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution coefficient (not used)
!>                 2) particle mass in AMU
!>                 3) safety factor: must be >1
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!>                 1) mpi rank
!> outputs:
!>   sup_pdf:      (real8) value of the probability density upper bound
function sup_pdf_current_density_uniform_phase(nx,x_min,x_max,fields,&
n_real_param,real_param,n_int_param,int_param) result(sup_pdf)
  use constants,          only: SPEED_OF_LIGHT,EL_CHG,MU_ZERO
  use mod_model_settings, only: n_var,var_zj
  use mod_fields,         only: fields_base
  implicit none
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: nx,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: sup_pdf
  !> Variables
  integer :: ifail
  real*8  :: density_tot,density_in,density_out
  real*8  :: pressure_tot,pressure_in,pressure_out
  real*8  :: kin_par_tot,kin_par_in,kin_par_out,mom_par_tot,mom_par_in
  real*8  :: mom_par_out
  real*8,dimension(2)     :: zj_minmax
  real*8,dimension(n_var) :: varmin,varmax
  real*8  :: max_pdf,min_pdf,sqrtpovermc2plus1_max,sqrtpovermc2plus1_min
  real*8  :: cos2pitch_max,cos2pitch_min
  !> Evalutate the upper extremum of the pdf
  call Integrals_3D(int_param(1),fields%node_list,fields%element_list,&
  density_tot,density_in,density_out,pressure_tot,pressure_in,pressure_out,&
  kin_par_tot,kin_par_in,kin_par_out,mom_par_tot,mom_par_in,mom_par_out,varmin,varmax)
  zj_minmax = [varmin(var_zj),varmax(var_zj)]
  call MPI_Bcast(zj_minmax,2,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  !> compute the maximum and the minimum of the pdf
  sqrtpovermc2plus1_max = sqrt((x_max(4)**2)/((real_param(2)*SPEED_OF_LIGHT)**2)+1.d0);
  cos2pitch_max = cos(x_max(5))**2;
  sqrtpovermc2plus1_min = sqrt((x_min(4)**2)/((real_param(2)*SPEED_OF_LIGHT)**2)+1.d0) 
  cos2pitch_min = cos(x_min(5))**2
  max_pdf = ((sqrtpovermc2plus1_max**3)-3.d0*sqrtpovermc2plus1_max) - &
            ((sqrtpovermc2plus1_min**3)-3.d0*sqrtpovermc2plus1_min);
  max_pdf = max_pdf*(cos2pitch_min - cos2pitch_max)*(x_max(6)-x_min(6))
  max_pdf =(-6.d0*real_param(3))/(max_pdf*x_min(7)*EL_CHG*MU_ZERO*(real_param(2)**3)*&
           (SPEED_OF_LIGHT**4)*x_min(1));
  min_pdf = max_pdf*zj_minmax(1); max_pdf = max_pdf*zj_minmax(2);
  !> check which between min_pdf and max_pdf has the maximum absolute value
  if(abs(max_pdf).ge.abs(min_pdf)) then
    sup_pdf = max_pdf
  else
    sup_pdf = min_pdf
  endif
end function sup_pdf_current_density_uniform_phase

!> Compute the gdf used for sampling the particle coordinates in 
!> phase space. The gdf used here is a uniform distribution in
!> cylindrical coordinates for the spatial coordinates and 
!> uniform spherical distribution for the momentum coordinates
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
function gdf_uniform_phase(nx,x,st,time,i_elm,fields,x_min,x_max,&
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
  !> Evalutate pdf
  gdf = 6.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))*&
  (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)*&
  (cos(x_min(5))-cos(x_max(5)))*(x_max(6)-x_min(6)))
  if(n_real_param.gt.0) gdf = real_param(1)*gdf
end function gdf_uniform_phase

!> Upper bound of the uniform phase space sampler distribution
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
!>   sup_gdf: (real8) value of the sampler probability density upper bound
function sup_gdf_uniform_phase(nx,x_min,x_max,n_real_param,real_param,&
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
  !> Evalutate the upper extremum of the pdf
  sup_gdf = 6.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))*&
  (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)*&
  (cos(x_min(5))-cos(x_max(5)))*(x_max(6)-x_min(6)))
end function sup_gdf_uniform_phase

!> GDF uniform sampler generating particle positions in phase space
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random numbers in [0,1)
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
!>   i_elm:        (integer) updated jorek mesh element number
!>   st:           (real8)(2) updated local mesh coordinates
!>   x:            (real8)(nx) particle position to accept
subroutine gdf_uniform_sampler(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Inputs-Outputs:
  integer,intent(inout)                       :: i_elm
  real*8,dimension(2),intent(inout)           :: st
  real*8,dimension(nx),intent(inout)          :: x
  !> Compute new particle position in phase space
  x(1) = sqrt(x_min(1)**2 + (x_max(1)**2 - x_min(1)**2)*x(1))
  x(2:3) = x_min(2:3) + (x_max(2:3)-x_min(2:3))*x(2:3)
  x(4) = (x_min(4)**3 + (x_max(4)**3-x_min(4)**3)*x(4))**(1d0/3d0)
  x(5) = acos(cos(x_min(5))-(cos(x_max(5))-cos(x_min(5)))*x(5))
  x(6:7) = x_min(6:7) + (x_max(6:7)-x_min(6:7))*x(6:7)
  !> find particle RZ coordinates
  call find_RZ(fields%node_list,fields%element_list,x(1),x(2),&
  x(1),x(2),i_elm,st(1),st(2),ifail)
end subroutine gdf_uniform_sampler

!> Dummy particle weight equal to 1 
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
!>                 1) plasma volume in SI units
!>                 2) total number of particles 
!>   n_int_param:  (integer) N# of integer input parameters
!>   int_param:    (integer)(n_int_param) integer weight parameters
!> outputs:
!>   weight:       (real) particle weight
function weight_uniform_one(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(weight)
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
  !> Compute particle weight
  weight = 1.d0
end function weight_uniform_one

!> Method used for computing the weight of each particle. Given that the pdf
!> is directly sampled via accept-reject method, the particle weight is 
!> considered uniform for all particles and equal to the total number of
!> physical particles (from the plasma volume) divided the total number of
!> simulated markers. A uniform distribution in the velocity space is used. 
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
!>                 1) plasma volume in SI units
!>                 2) total number of particles 
!>   n_int_param:  (integer) N# of integer input parameters
!>   int_param:    (integer)(n_int_param) integer weight parameters
!> outputs:
!>   weight:       (real) particle weight
function weight_uniform(nx,x,st,time,i_elm,fields,x_min,x_max,n_real_param,&
real_param,n_int_param,int_param) result(weight)
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
  !> Compute particle weight
  weight = (2.d0*real_param(1))/(real_param(2)*(x_max(1)**2-x_min(1)**2)*&
  (x_max(2)-x_min(2))*(x_max(3)-x_min(3)))
end function weight_uniform

!> Method used for computing the weight of each particle. Given that the pdf
!> is directly sampled via accept-reject method, the particle weight is 
!> considered uniform for all particles and equal to the total number of
!> physical particles (from the plasma current) divided the total number of
!> simulated markers. A uniform distribution in the velocity space is used. 
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
!>                 1) 3D integral of the plasma density in SI units
!>                 2) total number of particles
!>                 3) mass in AMU 
!>   n_int_param:  (integer) N# of integer input parameters
!>   int_param:    (integer)(n_int_param) integer weight parameters
!> outputs:
!>   weight:       (real) particle weight
function particle_weight_current_density_uniform_phase(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param) result(weight)
  use constants,  only: EL_CHG,SPEED_OF_LIGHT
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
  !> Variables:
  real*8 :: DUMMY_DOUBLE_1,DUMMY_DOUBLE_2
  !> Compute the particle weight
  DUMMY_DOUBLE_1 = sqrt((x_max(4)**2)/((real_param(3)*SPEED_OF_LIGHT)**2)+1.d0)
  DUMMY_DOUBLE_2 = sqrt((x_min(4)**2)/((real_param(3)*SPEED_OF_LIGHT)**2)+1.d0)
  weight = ((DUMMY_DOUBLE_1**3)-3.d0*DUMMY_DOUBLE_1) - ((DUMMY_DOUBLE_2**3)-3.d0*DUMMY_DOUBLE_2)
  weight = weight*(cos(x_min(5)+cos(x_max(5))))
  weight = (2.d0*real_param(1)*(x_max(4)**3 - x_min(4)**3))/&
  (weight*real_param(2)*x(7)*EL_CHG*(real_param(3)**3)*(SPEED_OF_LIGHT**4))
end function particle_weight_current_density_uniform_phase

!> Compute the total number of physical particles for the weight uniform pdf
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
!>                 1) plasma volume in SI units
!>                 2) total number of particles 
!>   n_int_param:  (integer) N# of integer input parameters
!>   int_param:    (integer)(n_int_param) integer weight parameters
!> outputs:
!>   n_phys_part:  (real8) number of physical particles
function n_physical_particle_weight_uniform(nx,time,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param) result(n_phys_part)
  use constants,  only: EL_CHG,SPEED_OF_LIGHT
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: nx,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: n_phys_part
  n_phys_part = real_param(1)
end function n_physical_particle_weight_uniform

!> Compute the total number of physical particles of the current density -
!> uniform momentum space distribution
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
!>                 1) 3D integral of the plasma density in SI units
!>                 2) total number of particles
!>                 3) mass in AMU 
!>   n_int_param:  (integer) N# of integer input parameters
!>   int_param:    (integer)(n_int_param) integer weight parameters
!> outputs:
!>   n_phys_part:  (real8) number of physical particles
function n_physical_particle_current_density_uniform_phase(nx,time,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param) result(n_phys_part)
  use constants,  only: EL_CHG,SPEED_OF_LIGHT
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: nx,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: n_phys_part
  !> Variables:
  real*8 :: DUMMY_DOUBLE_1,DUMMY_DOUBLE_2
  !> Compute the particle weight
  DUMMY_DOUBLE_1 = sqrt((x_max(4)**2)/((real_param(3)*SPEED_OF_LIGHT)**2)+1.d0)
  DUMMY_DOUBLE_2 = sqrt((x_min(4)**2)/((real_param(3)*SPEED_OF_LIGHT)**2)+1.d0)
  n_phys_part = ((DUMMY_DOUBLE_1**3)-3.d0*DUMMY_DOUBLE_1) - ((DUMMY_DOUBLE_2**3)-3.d0*DUMMY_DOUBLE_2)
  n_phys_part = n_phys_part*(cos(x_min(5)+cos(x_max(5))))
  n_phys_part = (2.d0*real_param(1)*(x_max(4)**3 - x_min(4)**3))/&
  (n_phys_part*x_min(7)*EL_CHG*(real_param(3)**3)*(SPEED_OF_LIGHT**4))
end function n_physical_particle_current_density_uniform_phase

end program test_initialisation_phase_space
