!> Module containing the base functionality for the initialization of particles in 
!> configuration space (6D) by rejection sampling
module initialisers_base
  use mod_rng
  use data_structure
  use mod_particle_types
  use constants
  use mod_interp
  use mod_rej_f
  implicit none
  private
  public initialise_particles, adjust_particle_weights
  public set_velocity_from_T, domain_bounding_box, initialise_particles_H_mu_psi
  public initialise_particles_H_mu_psi_phiplanes
  public set_particle_weights_canonical_maxwellian, normalize_with_projection
  public weigh_with_interp_f
  public normalize_with_projection_at_gc
  public real_f,real_arr_inout_s,part_inout_s,initialise_particles_in_phase_space

  interface
    subroutine find_RZ(node_list,element_list,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
      use data_structure
      type (type_node_list), intent(in)    :: node_list
      type (type_element_list), intent(in) :: element_list
      real*8, intent(in)     :: R_find, Z_find
      real*8, intent(out)    :: R_out,Z_out,s_out,t_out
      integer, intent(inout) :: ielm_out
      integer, intent(out)   :: ifail
    end subroutine find_RZ
    function real_f(n_x,x,st,time,i_elm,fields,x_min,x_max,&
    n_real_param,real_param,n_int_param,int_param)
      use mod_fields, only: fields_base
      implicit none
      !> inputs:
      integer,intent(in)                                  :: n_x,i_elm
      integer,intent(in)                                  :: n_real_param,n_int_param
      integer,dimension(:),allocatable,intent(in)         :: int_param
      real*8,intent(in)                                   :: time
      real*8,dimension(n_x),intent(in)                    :: x,x_min,x_max
      real*8,dimension(2),intent(in)                      :: st
      real*8,dimension(:),allocatable,intent(in)          :: real_param
      class(fields_base),intent(in)                       :: fields
      !> outputs:
      real*8                                              :: real_f
    end function real_f
    subroutine real_arr_inout_s(n_x,x,st,time,i_elm,fields,x_min,x_max,&
    n_real_param,real_param,n_int_param,int_param)
      use mod_fields, only: fields_base
      implicit none
      !> inputs:
      integer,intent(in)                                  :: n_x
      integer,intent(in)                                  :: n_real_param,n_int_param
      integer,dimension(:),allocatable,intent(in)         :: int_param
      real*8,intent(in)                                   :: time
      real*8,dimension(n_x),intent(in)                    :: x_min,x_max
      real*8,dimension(:),allocatable,intent(in)          :: real_param
      class(fields_base),intent(in)                       :: fields
      !> inputs-outputs:
      integer,intent(inout)                               :: i_elm
      real*8,dimension(2),intent(inout)                   :: st
      real*8,dimension(n_x),intent(inout)                 :: x
    end subroutine real_arr_inout_s
    subroutine part_inout_s(p_inout,n_x,x,time,fields,n_real_param,&
    real_param,n_int_param,int_param)
      use mod_particle_types, only: particle_base
      use mod_fields,         only: fields_base
      implicit none
      !> inputs:
      integer,intent(in)               :: n_x
      integer,intent(in)               :: n_int_param,n_real_param
      integer,dimension(:),allocatable,intent(in) :: int_param
      real*8,intent(in)                :: time
      real*8,dimension(n_x),intent(in) :: x
      real*8,dimension(:),allocatable,intent(in)  :: real_param
      !> inputs-outputs:
      class(particle_base),intent(inout) :: p_inout
      class(fields_base),intent(in)      :: fields
    end subroutine
  end interface

contains
!> Set positions for particles by rejection sampling from geometric and mhd
!> variables after collecting with transform, within Rbound, Zbound and Phibound
!> if present. See [[test_rejection_sampling]] for examples.
subroutine initialise_particles(particles, node_list, element_list, &
  rng, space_pdf, Rbound, Zbound, Phibound, &
  rng_n_streams_round_off_in)
  use mpi
  use mod_sampling
  use mod_random_seed
  use mod_interp
!$ use omp_lib
  implicit none

  class(particle_base), dimension(:), intent(inout) :: particles
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  class(type_rng), intent(in)                       :: rng       !< What type of random number generator to use. Is re-seeded in the subroutine.
  type(spatial_pdf), intent(in), optional           :: space_pdf
  real*8, dimension(2), intent(in), optional        :: Rbound, Zbound, Phibound   !< Between which coordinates to sample (RZPhi).
                                                                                  !< if omitted, determine automatically from node_list
  logical, intent(in), optional                     :: rng_n_streams_round_off_in !< round-off the rng n_streams at 2**ceil

  ! Internal variables
  real*8  :: R, Z, phi, s, t, DUMMY_REAL
  real*8  :: Rbox(2), Zbox(2), Phibox(2)
  integer :: i, j, ifail
  real*8  :: ran(7)
  integer :: i_elm
  real*8  :: t0, t1, ostart, oend
  integer :: seq, n_streams, n_threads, i_thread
  integer :: my_id, n_mpi
  integer :: seed
  logical :: rng_n_streams_round_off, use_rej
  class(type_rng), allocatable, dimension(:) :: rngs ! The RNGs for all the threads
  integer, dimension(:), allocatable         :: i_to_find
  logical, dimension(:), allocatable         :: not_found

  ostart = 0.d0
  oend   = 0.d0

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ifail)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ifail)
  rng_n_streams_round_off = .false.
  if(present(rng_n_streams_round_off_in)) rng_n_streams_round_off = rng_n_streams_round_off_in

  use_rej = .false.
  if (present(space_pdf)) then
    call check_spatial_pdf(space_pdf, "initialise_particles")
    use_rej = associated(space_pdf%f)
  endif

  ! Setup bounding boxes
  call domain_bounding_box(node_list, element_list, Rbox(1), Rbox(2), Zbox(1), Zbox(2))
  Phibox = [0.d0, TWOPI]
  ! Check if the requested bounding boxes have any overlap with the domain
  ! Does not check for combinations of R and Z
  if (present(Rbound)) then
    if (maxval(Rbound) .gt. minval(Rbox) .and. maxval(Rbox) .gt. minval(Rbound)) then
      Rbox = Rbound
    else
      write(*,*) "ERROR: no overlap between domain and requested bounding box in R, domain=", Rbox, ", box=", Rbound
      write(*,*) "Sampling from whole domain in R"
    end if
  end if
  if (present(Zbound)) then
    if (maxval(Zbound) .gt. minval(Zbox) .and. maxval(Zbox) .gt. minval(Zbound)) then
      Zbox = Zbound
    else
      write(*,*) "ERROR: no overlap between domain and requested bounding box in Z, domain=", Zbox, ", box=", Zbound
      write(*,*) "Sampling from whole domain in Z"
    end if
  end if
  if (present(Phibound)) PhiBox = Phibound

  ! Calculate a single random seed and communicate it over MPI
  if (my_id .eq. 0) seed = random_seed()
  call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ifail)

  ! Prepare list of particles to seed
  allocate(i_to_find(size(particles,1)),not_found(size(particles,1)))
  i_to_find = [(i, i=1,size(particles,1))] ! which particles still to do
  not_found = .true. ! whether this one has been sampled succesfully

  ! Setup (Q)RNGs, one per thread
  n_threads = 1
!$ n_threads = omp_get_max_threads()
  allocate(rngs(0:n_threads-1), source=rng)
  n_streams = n_mpi*n_threads ! Works only for homogeneous environments!

  do i_thread=0,n_threads-1
    seq = my_id*n_threads + i_thread + 1
    call rngs(i_thread)%initialize(7, seed, n_streams, seq, &
    ifail, round_off_n_streams_in=rng_n_streams_round_off)
    if (ifail .ne. 0) call MPI_ABORT(MPI_COMM_WORLD, -1, ifail)
  end do

  call cpu_time(t0)
!$ ostart = omp_get_wtime()

  ! Filter over all particles to sample them, repeat for rejected positions until
  ! empty. This is required if the distribution has some correlation with the
  ! samples mod something (as is the case for the Sobol sequence).
  ! Even then, an inbalance in openmp scheduling or number of particles per node
  ! could cause a slight correlation
  ! in the output. This could be worse if the distribution is very narrow.
  ! In that case the Sobol series should also be implemented in 64-bits as the total
  ! number of values is 2^31 now.
  ! TODO fix also for MPI or broadcast to nodes
  do while (any(not_found))
    ! default(shared) is very dangerous but needed due to gfortran failures.
    ! be very careful (error message for default(none) below)
    ! Error: ‘__vtab_mod_particle_types_Particle_kinetic_leapfrog’ not specified in enclosing ‘parallel’
#ifdef __GFORTRAN__
    !$omp parallel default(shared) &
#else
    !$omp parallel default(none) &
#endif
    !$omp   shared(particles, node_list, element_list, Rbox, Zbox, PhiBox, space_pdf, use_rej, &
    !$omp          rngs, n_threads, n_streams, seed, my_id, i_to_find, not_found) &
    !$omp   private(j, i, R, Z, phi, i_elm, s, t, ifail, seq, ran, i_thread, DUMMY_REAL)
    i_thread = 0
!$  i_thread=omp_get_thread_num()
    !$omp do schedule(static)
    do i=1,size(i_to_find,1)
      j = i_to_find(i)
      ! Generate a random position to put this particle
      call rngs(i_thread)%next(ran)
      call transform_uniform_cylindrical(ran(1:3), Rbox, Zbox, PhiBox, R, Z, phi)

      call find_RZ(node_list,element_list,R,Z,DUMMY_REAL,DUMMY_REAL,i_elm,s,t,ifail)
      if (ifail .eq. 0) then
      !> Reject this position if it fails the spatial pdf
        if (use_rej) then
          if (ran(4) .ge. eval_rej_f(node_list, element_list, i_elm, s, t, phi, R, Z, space_pdf)) cycle
        end if
        particles(j)%x     = [R, z, phi]
        particles(j)%i_elm = i_elm
        particles(j)%st    = [s, t]
        select type(pa => particles(j))
          type is (particle_kinetic_leapfrog)
            pa%v = ran(5:7) !> Save other components of this point for velocity init in a later routine
        endselect
        not_found(i) = .false.
      end if
    enddo
    !$omp end do
    !$omp end parallel
    ! now pack only the indices of particles we still need to do
    i_to_find = pack(i_to_find, not_found) ! implicitly allocates
    deallocate(not_found); allocate(not_found(size(i_to_find,1)))
    not_found = .true.
  end do

  call cpu_time(t1)
!$ oend = omp_get_wtime()
  write(*,'(i5,A,2f12.4)') my_id, ' Time particle initialize cpu/wall :',t1-t0, oend-ostart
  if (my_id .eq. 0) then
    write(*,*) '* done initialising particles    *'
    write(*,*) '**********************************'
  endif
end subroutine initialise_particles

!> Initialise particle in phase space given a generic distribution in space and momentum
!> space. The acceptance-rejection method is used for generating the particle population.
!> The particle population is sampled from a generic distribution probability function
!> gdf. Both the gdf and the gdf sampling routines must be provided as inputs.
!> Check examples/test_initialisation_phase_space.f90 for an example of its usage.
!> Inputs:
!>   n_variables:                (integer) dimensionality of the phase space
!>   particles:                  (particle_base)(n_particles) particle array to be initialised
!>   fields:                     (fields_base) jorek MHD fields data structure
!>   rng_base:                   (type_rng) type of random number generator to be used
!>   pdf:                        (real_f) particle probability density function
!>   weight_f:                   (real_f) method computing the particle weight
!>   gdf:                        (real_f) probability density function used for sampling 
!>                               the particle coordinates
!>   gdf_sampler:                (real_arr_inout_s) method for generating samples of
!>                               particle coordinates from the gdf probability density
!>   sup_pdf:                    (real8) pdf upper extremum
!>   sup_gdf:                    (real8) gdf upper extremum
!>   sample_to_particle:         (part_inout_s) method use for transforming a sample in
!>                               sample coordinates in particle coordinates
!>   mass:                       (real8) particle mass
!>   time:                       (real8) physical time at which particles are initialised
!>   phase_bounds:               (real8)(n_variables,2) minima (first column) and maxima 
!>                               (second column) of the sampling phase space interval
!>   n_real_pdf_param_in:        (integer) number of real parameters of the pdf
!>   real_pdf_param_in:          (real8)(n_real_pdf_param_in) pdf real parameters
!>   n_int_pdf_param_in:         (integer) number of integer parameters of the pdf
!>   int_pdf_param_in:           (integer)(n_int_pdf_param_in) pdf integer parameters
!>   n_real_weight_param_in:     (integer) number of real parameters of the weight
!>   real_weight_param_in:       (real8)(n_real_weight_param_in) weight real parameters
!>   n_int_weight_param_in:      (integer) number of integer parameters of the weight
!>   int_weight_param_in:        (integer)(n_int_weight_param_in) weight integer parameters
!>   n_real_gdf_param_in:        (integer) number of real parameters of the gdf
!>   real_gdf_param_in:          (real8)(n_real_gdf_param_in) gdf real parameters
!>   n_int_gdf_param_in:         (integer) number of integer parameters of the gdf
!>   int_gdf_param_in:           (integer)(n_int_gdf_param_in) gdf integer parameters
!>   n_real_samp_to_part_in:     (real8) size of the real parameter array of the
!>                               sample to particle method
!>   real_samp_to_part_in:       (real8) intger parameter array od the sample to particle method
!>   n_int_samp_to_part_in:      (integer) size of the integer parameter array of the
!>                               sample to particle method
!>   int_samp_to_part_in:        (integer) intger parameter array od the sample to particle method
!>   rng_n_streams_round_off_in: (logical) round-off the rng n_streams at 2**ceil

!> Outputs:
!>   particles:              (particle_base)(n_particles) initialised particle array
subroutine initialise_particles_in_phase_space(n_variables, particles, fields, rng_base, pdf, &
  weight_f, gdf, gdf_sampler, sup_pdf, sup_gdf, sample_to_particle, mass, time, phase_bounds, &
  n_real_pdf_param_in, real_pdf_param_in, n_int_pdf_param_in, int_pdf_param_in, &
  n_real_weight_param_in, real_weight_param_in, n_int_weight_param_in, int_weight_param_in, &
  n_real_gdf_param_in, real_gdf_param_in, n_int_gdf_param_in, int_gdf_param_in, &
  n_real_samp_to_part_param_in,real_samp_to_part_param_in,n_int_samp_to_part_param_in, &
  int_samp_to_part_param_in,rng_n_streams_round_off_in)
  use mod_fields,                only: fields_base
  use mod_random_seed,           only: random_seed
  use mod_particle_types,        only: particle_base
  use mod_rng
!$ use omp_lib
  
  !> parameters
  integer,parameter :: chunksize=64

  !> inputs-outputs
  class(particle_base), dimension(:), intent(inout) :: particles
  !> inputs
  class(fields_base),intent(in)                        :: fields
  class(type_rng),intent(in)                           :: rng_base !< What type of random number generator to use (will be reseeded here)
  procedure(real_f)                                    :: pdf,weight_f,gdf
  procedure(real_arr_inout_s)                          :: gdf_sampler
  procedure(part_inout_s)                              :: sample_to_particle
  integer,intent(in)                                   :: n_variables
  integer,intent(in),optional                          :: n_real_pdf_param_in,n_int_pdf_param_in
  integer,intent(in),optional                          :: n_real_weight_param_in,n_int_weight_param_in
  integer,intent(in),optional                          :: n_real_gdf_param_in,n_int_gdf_param_in
  integer,intent(in),optional                          :: n_real_samp_to_part_param_in
  integer,intent(in),optional                          :: n_int_samp_to_part_param_in 
  integer,dimension(:),allocatable,intent(in),optional :: int_pdf_param_in,int_weight_param_in
  integer,dimension(:),allocatable,intent(in),optional :: int_gdf_param_in,int_samp_to_part_param_in
  real*8,intent(in)                                    :: mass,time,sup_pdf,sup_gdf
  real*8,dimension(n_variables,2),intent(in)           :: phase_bounds
  real*8,dimension(:),allocatable,intent(in),optional  :: real_pdf_param_in,real_weight_param_in
  real*8,dimension(:),allocatable,intent(in),optional  :: real_gdf_param_in,real_samp_to_part_param_in
  logical,                        intent(in),optional  :: rng_n_streams_round_off_in !< round-off the rng n_streams at 2**ceil

  !> internal variables
  class(type_rng),dimension(:),allocatable :: rngs 
  integer                                  :: t0,t1,my_id,n_mpi,n_threads,thread_id,ifail
  integer                                  :: ii,jj,n_particles,i_elm
  integer                                  :: n_real_pdf_param,n_int_pdf_param
  integer                                  :: n_real_weight_param,n_int_weight_param
  integer                                  :: n_real_gdf_param,n_int_gdf_param
  integer                                  :: n_real_samp_to_part_param,n_int_samp_to_part_param
  integer,dimension(:),allocatable         :: int_pdf_param,int_weight_param,int_gdf_param
  integer,dimension(:),allocatable         :: int_samp_to_part_param
  logical                                  :: rng_n_streams_round_off
  real*8                                   :: weight,one_over_sup_pdf,one_over_sup_gdf
  real*8,dimension(2)                      :: st
  real*8,dimension(:),allocatable          :: variables,real_pdf_param,real_weight_param
  real*8,dimension(:),allocatable          :: real_gdf_param,real_samp_to_part_param

  !> extract id and size of the MPI Communicator
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ifail)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ifail)
  rng_n_streams_round_off = .false.
  if(present(rng_n_streams_round_off_in)) rng_n_streams_round_off = rng_n_streams_round_off_in
  !> initialize random number generator
  n_threads = 1
!$ n_threads = omp_get_max_threads()
  allocate(variables(n_variables+1))
  allocate(rngs(n_threads),source=rng_base)
  do ii=1,n_threads
    call rngs(ii)%initialize(n_variables+1, random_seed(), n_mpi*n_threads, &
    my_id*n_threads+ii,ifail,round_off_n_streams_in=rng_n_streams_round_off)
    if (ifail.ne.0) call MPI_ABORT(MPI_COMM_WORLD, -1, ifail)
  end do

  !> Initialise variables needed in the loop
  n_particles = size(particles);
  one_over_sup_pdf = 1d0/sup_pdf; one_over_sup_gdf = 1d0/sup_gdf;
  n_real_pdf_param = 0; if(present(n_real_pdf_param_in)) n_real_pdf_param = n_real_pdf_param_in
  if((present(real_pdf_param_in)).and.(n_real_pdf_param.gt.0)) then
    if(allocated(real_pdf_param_in)) then
      allocate(real_pdf_param(n_real_pdf_param)); real_pdf_param = real_pdf_param_in;
    endif
  endif
  n_int_pdf_param = 0; if(present(n_int_pdf_param_in)) n_int_pdf_param = n_int_pdf_param_in;
  if((present(int_pdf_param_in)).and.(n_int_pdf_param.gt.0)) then
    if(allocated(int_pdf_param_in)) then
      allocate(int_pdf_param(n_int_pdf_param)); int_pdf_param = int_pdf_param_in;
    endif
  endif
  n_real_weight_param = 0; if(present(n_real_weight_param_in)) n_real_weight_param = n_real_weight_param_in
  if((present(real_weight_param_in)).and.(n_real_weight_param.gt.0)) then
    if(allocated(real_weight_param_in)) then
      allocate(real_weight_param(n_real_weight_param)); real_weight_param = real_weight_param_in;
    endif
  endif
  n_int_weight_param = 0; if(present(n_int_weight_param_in)) n_int_weight_param = n_int_weight_param_in;
  if((present(int_weight_param_in)).and.(n_int_weight_param.gt.0)) then
    if(allocated(int_weight_param_in)) then
      allocate(int_weight_param(n_int_weight_param)); int_weight_param = int_weight_param_in;
    endif
  endif
  n_real_gdf_param = 0; if(present(n_real_gdf_param_in)) n_real_gdf_param = n_real_gdf_param_in
  if((present(real_gdf_param_in)).and.(n_real_gdf_param.gt.0)) then
    if(allocated(real_gdf_param_in)) then
      allocate(real_gdf_param(n_real_gdf_param)); real_gdf_param = real_gdf_param_in;
    endif
  endif
  n_int_gdf_param = 0; if(present(n_int_gdf_param_in)) n_int_gdf_param = n_int_gdf_param_in;
  if((present(int_gdf_param_in)).and.(n_int_gdf_param.gt.0)) then
    if(allocated(int_gdf_param_in)) then
      allocate(int_gdf_param(n_int_gdf_param)); int_gdf_param = int_gdf_param_in;
    endif
  endif
  n_real_samp_to_part_param = 0; if(present(n_real_samp_to_part_param_in)) &
  n_real_samp_to_part_param = n_real_samp_to_part_param_in
  if((present(real_samp_to_part_param_in)).and.(n_real_samp_to_part_param.gt.0)) then
    if(allocated(real_samp_to_part_param_in)) then
      allocate(real_samp_to_part_param(n_real_samp_to_part_param)); 
      real_samp_to_part_param = real_samp_to_part_param_in
    endif
  endif
  n_int_samp_to_part_param = 0; if(present(n_int_samp_to_part_param_in)) &
  n_int_samp_to_part_param = n_int_samp_to_part_param_in;
  if((present(int_samp_to_part_param_in)).and.(n_int_samp_to_part_param.gt.0)) then
    if(allocated(int_samp_to_part_param)) then
      allocate(int_samp_to_part_param(n_int_samp_to_part_param)); int_samp_to_part_param = int_samp_to_part_param_in;
    endif
  endif
  call system_clock(t0)

  !> Loop on the particles
#ifndef __NVCOMPILER
    !$omp parallel default(shared) &
    !$omp firstprivate(n_variables,n_particles,mass,time,phase_bounds,one_over_sup_pdf,&
    !$omp one_over_sup_gdf,n_real_pdf_param,n_int_pdf_param,real_pdf_param,&
    !$omp int_pdf_param,n_real_weight_param,n_int_weight_param,real_weight_param,&
    !$omp int_weight_param,n_real_gdf_param,n_int_gdf_param,real_gdf_param,int_gdf_param,&
    !$omp n_real_samp_to_part_param,n_int_samp_to_part_param,real_samp_to_part_param,&
    !$omp int_samp_to_part_param) &
    !$omp private(ii,variables,thread_id,i_elm,st,weight,ifail)
    thread_id = 1
    !$ thread_id = omp_get_thread_num()+1
    !$omp do schedule(dynamic,chunksize)
#endif
    do ii=1,n_particles
      i_elm = 0; st = -1.d0; weight = 0d0;
      !> loop until the particle is not valid, it can slow down the code
      !> but before trying a manual load balacing has done in initialise_particles_H_mu_psi
      !> let's check how the openMP dynamic scheduling performs using different chunksize
      do while(rejection_funct_gpdf(n_variables,variables(1:n_variables),st,time,i_elm,weight,&
        variables(n_variables+1),phase_bounds(:,1),phase_bounds(:,2),fields,pdf,gdf,&
        one_over_sup_pdf,one_over_sup_gdf,n_real_pdf_param,real_pdf_param,n_int_pdf_param,&
        int_pdf_param,n_real_gdf_param,real_gdf_param,n_int_gdf_param,int_gdf_param))
        !> uniform sampling in cylindrical coordinates for the physical space (R,Z,phi),
        !> and general coordinates for the phase space
        call rngs(thread_id)%next(variables)
        call gdf_sampler(n_variables,variables(1:n_variables),st,time,i_elm,fields,&
        phase_bounds(:,1),phase_bounds(:,2),n_real_gdf_param,real_gdf_param,&
        n_int_gdf_param,int_gdf_param)
        weight = weight_f(n_variables,variables(1:n_variables),st,time,&
                 i_elm,fields,phase_bounds(:,1),phase_bounds(:,2),n_real_weight_param,&
                 real_weight_param,n_int_weight_param,int_weight_param)
      enddo
      !> store the values of accepted particles 
      particles(ii)%x      = variables(1:3)
      particles(ii)%st     = st
      particles(ii)%i_elm  = i_elm
      particles(ii)%weight = weight
      !> transform the remaining variables of the sampe in particle coordinates
      call sample_to_particle(particles(ii),n_variables,variables(1:n_variables),time,fields,&
      n_real_samp_to_part_param,real_samp_to_part_param,n_int_samp_to_part_param,&
      int_samp_to_part_param)
    enddo
#ifndef __NVCOMPILER
    !$omp end do
    !$omp end parallel
#endif

  !> clean-up
  call system_clock(t1)
  deallocate(variables); deallocate(rngs);
  if(allocated(real_pdf_param))          deallocate(real_pdf_param)
  if(allocated(int_pdf_param))           deallocate(int_pdf_param)
  if(allocated(real_weight_param))       deallocate(real_weight_param)
  if(allocated(int_weight_param))        deallocate(int_weight_param)
  if(allocated(real_gdf_param))          deallocate(real_gdf_param)
  if(allocated(int_gdf_param))           deallocate(int_gdf_param)
  if(allocated(real_samp_to_part_param)) deallocate(real_samp_to_part_param)
  if(allocated(int_samp_to_part_param))  deallocate(int_samp_to_part_param)
  write(*,'(i5,A,2f12.4)') my_id, ' Time particle initialize system (s) :',real(t1-t0,kind=8)/1d3
  if (my_id .eq. 0) then
    write(*,*) '* done initialising particles    *'
    write(*,*) '**********************************'
  endif

end subroutine initialise_particles_in_phase_space

!> rejection function of the acceptance - rejection method
!> inputs:
!>   n_x:              (integer) number of variables
!>   x:                (real8)(n_x) variables
!>   st:               (real8)(2) jorek mesh local coordinates
!>   i_elm:            (integer) jorek element number
!>   weight:           (real8) particle weight
!>   rand:             (real8) uniformly distributed random number [0,1]
!>   fields:           (fields_base) jorek MHD fields
!>   pdf:              (real_f) procedure returning the value of the
!>                     probability density function at a given point
!>   gdf:              (real_f) procedure returning the value of the
!>                     probability density function used for generating
!>                     the variables to be tested
!>   one_over_sup_pdf: (real8) 1/upper bound of the pdf
!>   one_over_sup_gdf: (real8) 1/upper bound of the gdf
!>   n_real_pdf_param: (integer) N# of double parameters of the pdf
!>   real_pdf_param:   (real8)(n_real_pdf_param) pdf double parameters
!>   n_int_pdf_param:  (integer) N# of integer parameters of the pdf
!>   int_pdf_param:    (integer)(n_int_pdf_param) pdf integer parameters
!>   n_real_gdf_param: (integer) N# of double parameters of the gdf
!>   real_gdf_param:   (real8)(n_real_pdf_param) gdf double parameters
!>   n_int_gdf_param:  (integer) N# of integer parameters of the gdf
!>   int_gdf_param:    (integer)(n_int_pdf_param) gdf integer parameters
!> outputs:
!>   rej: (logical) if true the sample is rejected
function rejection_funct_gpdf(n_x,x,st,time,i_elm,weight,rand,&
x_min,x_max,fields,pdf,gdf,one_over_sup_pdf,one_over_sup_gdf,&
n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param,&
n_real_gdf_param,real_gdf_param,n_int_gdf_param,int_gdf_param) result(rej)
  use mod_fields, only: fields_base
  implicit none
  !> input variables
  class(fields_base),intent(in)    :: fields
  integer,intent(in)               :: n_x,i_elm
  integer,intent(in)               :: n_real_pdf_param, n_int_pdf_param
  integer,intent(in)               :: n_real_gdf_param, n_int_gdf_param
  integer,dimension(:),allocatable,intent(in) :: int_pdf_param,int_gdf_param
  real*8,intent(in)                :: weight,time,rand,one_over_sup_pdf,one_over_sup_gdf
  real*8,dimension(2)              :: st
  real*8,dimension(n_x),intent(in) :: x,x_min,x_max
  real*8,dimension(:),allocatable,intent(in) :: real_pdf_param,real_gdf_param
  procedure(real_f)                :: pdf,gdf
  !> output variables
  logical :: rej
  real*8 :: norm_pdf,norm_gdf

  !> check if the particle is valid
  rej = .true.
  if(i_elm.le.0) return
  if(weight.le.0d0) return
  if((st(1).lt.0.d0).or.(st(1).gt.1.d0)) return
  if((st(2).lt.0.d0).or.(st(2).gt.1.d0)) return
  !> check if the pdf is valid
  norm_pdf = one_over_sup_pdf*pdf(n_x,x,st,time,i_elm,fields,x_min,&
  x_max,n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
  norm_gdf = one_over_sup_gdf*gdf(n_x,x,st,time,i_elm,fields,x_min,&
  x_max,n_real_gdf_param,real_gdf_param,n_int_gdf_param,int_gdf_param)
  if(norm_pdf.lt.0d0) then
    write(*,*) 'pdf/sup(pdf) smaller than 0: skip!' 
    return
  endif
  if(norm_gdf.lt.0d0) then
    write(*,*) 'gdf/sup(gdf) smaller than 0: skip!' 
    return
  endif
  if(norm_pdf.gt.1d0) write(*,*) 'WARNING: normalised pdf > 1: increase the pdf extremum safety factor'
  if(norm_gdf.gt.1d0) write(*,*) 'WARNING: normalised gdf > 1: increase the pdf extremum safety factor'
  !> reject or accept solution
  if((rand*norm_gdf).le.norm_pdf) rej = .false.
  
end function rejection_funct_gpdf

!> Initialise particle positions in E, mu, (psi, theta|R, Z), phi, gamma (gyrophase) space.
!> Set Psi_transform to transform from [0,1] to your desired range
subroutine initialise_particles_H_mu_psi(particles, fields, rng_base, mass, T_maxwell, &
  Theta_transform, Psi_transform, alpha, E_max, include_vpar, uniform_space, &
  space_pdf, cor, charge, rng_n_streams_round_off_in)
  use mod_rng
  use mod_fields
  use mod_random_seed
  use mod_sampling
  use constants
  use phys_module, only: F0, central_density
  use mod_coronal
  use mod_boris, only: gc_to_kinetic_leapfrog, kinetic_to_kinetic_leapfrog, gc_to_kinetic
  use mod_gc_variational, only: convert_gc_to_gc_vpar
  use mpi
  use mod_interp
  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  class(fields_base),    intent(in)                 :: fields
  class(type_rng),       intent(in)                 :: rng_base !< What type of random number generator to use (will be reseeded here)
  real*8,                intent(in)                 :: mass
  real*8,                external,   optional       :: Theta_transform !< Function to transform 0-1 to the theta-domain
  real*8,                external,   optional       :: Psi_transform !< Function to transform 0-1 to the Psi-domain
  !< if omitted, determine automatically from node_list
  real*8,                intent(in), optional       :: alpha !< Make more fast (>0) or slow (<0) particles and weigh them appropriately
  real*8,                intent(in), optional       :: E_max !< If alpha=1 we select particles from a block-distribution, up to E_max
  logical,               intent(in), optional       :: include_vpar !< Initialize particles with local parallel velocity
  logical,               intent(in), optional       :: uniform_space !< Do not
  !< use {psi,theta}_transform if present bute rejection sampling in RZ
  type(spatial_pdf),     intent(in), optional       :: space_pdf !< spatial pdf to reject against when uniform_space is set. If absent the RZ sampling is uniform
  type(coronal),         intent(in), optional       :: cor !< Coronal equilibrium datatype for this particle. If unset, do not alter q
  integer,               intent(in), optional       :: charge !< Use this if cor is not present
  real*8,                intent(in), optional       :: T_Maxwell !< constant Maxwellian temperature [eV]
  logical,               intent(in), optional       :: rng_n_streams_round_off_in !< round-off the rng n_streams at 2**ceil

  ! Internal variables
  type(particle_gc)      :: particle
  type(particle_kinetic) :: particle_kinetic_tmp
  class(type_rng), allocatable :: rng
  real*8  :: ran(8)
  real*8  :: H, muB, chi, V2, v_par
  real*8  :: psi, psimin, psimax, theta, phi
  real*8  :: R, Z, inv_st_jac, psi_r, psi_z, B(3)
#ifdef fullmhd
  real*8    :: A3, AR, AZ, A3_R, A3_Z, AR_Z, AR_p, AZ_R, AZ_P, Fprof
  real*8, dimension(3)                :: P, P_s, P_t, P_phi
#else 
  real*8, dimension(1)                :: P, P_s, P_t, P_phi
#endif
  
  real*8  :: R_s, R_t, Z_s, Z_t
  real*8  :: s, t, u_init_max, temp, u
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
  integer :: i_elm, i, j, ifail, my_id, n_mpi, ierr
  real*8, dimension(fields%element_list%n_elements,2)    :: psi_minmax_list
  real*8, allocatable, dimension(:,:)             :: rans
  class(particle_base), dimension(:), allocatable :: particles_tmp
  logical, dimension(:), allocatable              :: found
  real*8  :: t0, t1
  real*8  :: Rbox(2), Zbox(2), DUMMY_R, DUMMY_Z
  integer :: blocksize, prev_blocksize, particles_to_do_local, particles_done_local
  integer :: to_find, n_tries_now, n_found
  logical :: all_done, init_uniform_space, my_include_vpar, rng_n_streams_round_off, use_rej
  real*8  :: my_alpha

  rng_n_streams_round_off = .false.
  if(present(rng_n_streams_round_off_in)) rng_n_streams_round_off = rng_n_streams_round_off_in

  init_uniform_space = .false.
  if (present(uniform_space) .and. uniform_space) then
    init_uniform_space = .true.
    call domain_bounding_box(fields%node_list, fields%element_list, Rbox(1), Rbox(2), Zbox(1), Zbox(2))
  end if

  my_include_vpar = .false.
  if (present(include_vpar)) my_include_vpar = include_vpar
  ! Check if we are in the 3,4,5 series of models
  if (my_include_vpar .and. .not. with_Vpar) then
    write(*,*) "WARNING: This model, ", jorek_model, "does not support parallel flows, disabling"
    my_include_vpar = .false.
  end if

  use_rej = .false.
  if (present(space_pdf)) then
    call check_spatial_pdf(space_pdf, "initialise_particles_H_mu_psi")
    use_rej = associated(space_pdf%f)
  endif

  if (present(alpha)) then
    write(*,*) "alpha not implemented yet"
    call exit(1)
    if (alpha .eq. 1.d0 .and. .not. present(E_max)) then
      write(*,*) "E_max is required if alpha=1"
      call exit(1)
    end if
    my_alpha = alpha
  else
    my_alpha = 0.d0
  end if


  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ierr)

  if (.not. init_uniform_space) then
    psimin= 1d10
    psimax=-1d10
    ! Preparatory work: determine psi_min,max
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none)        &
    !$omp shared(fields,  psi_minmax_list) &
#endif
    !$omp private(i_elm) reduction(min:psimin) &
    !$omp reduction(max:psimax)
    do i_elm=1, fields%element_list%n_elements
      call psi_minmax(fields%node_list, fields%element_list, i_elm, psi_minmax_list(i_elm,1), psi_minmax_list(i_elm,2))
      psimin = min(psi_minmax_list(i_elm,1),psimin)
      psimax = max(psi_minmax_list(i_elm,2),psimax)
    end do
    !$omp end parallel do
  end if

  ! Preparatory work: setup RNG
  allocate(rng,source=rng_base)
  call rng%initialize(8, random_seed(), 1, 1, ifail, round_off_n_streams_in=rng_n_streams_round_off)

  ! Preparatory work: get R_axis, Z_axis
  if (.not. init_uniform_space) then
    call find_axis(my_id, fields%node_list, fields%element_list, psi_axis, R_axis, Z_axis, i_elm, s_axis, t_axis, ifail)
  end if

  ! Preset i_elm to 0 so that by default particles are lost
  particles(:)%i_elm = 0

  ! Assume equal number of particles everywhere
  particles_to_do_local  = size(particles)
  particles_done_local   = 0
  prev_blocksize = 0 ! initial block size
  all_done = .false.

  do while (.not. all_done)
    to_find = (particles_to_do_local-particles_done_local)
    n_tries_now = to_find
    ! Add a little bit extra to cut off tail of 1/x^n
    if (n_tries_now .lt. 64 .and. n_tries_now .gt. 0) n_tries_now = 64

    ! Calculate starting-point for this block
    ! random numbers are generated in blocks per mpi process
    ! blocks are grouped in a superblock of block*n_mpi size
    ! blocks must be equal size
    ! we communicate a blocksize here, and keep a running index of the previous starting points
    ! make this the biggest of n_tries_now
    call MPI_AllReduce(n_tries_now, blocksize, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    call rng%jump_ahead((n_mpi-my_id-1)*prev_blocksize + my_id*blocksize)

    prev_blocksize = blocksize

    ! Generate the random numbers
    call cpu_time(t0)
    allocate(rans(8,blocksize))
    do i=1,blocksize
      call rng%next(rans(:,i))
    end do
    call cpu_time(t1)

    ! Allocate some temporary storage
    select type (particles)
      type is (particle_kinetic)
        write(*,*) "ERROR: particle_kinetic not supported yet for initialize_particles_H_mu_psi"
      type is (particle_kinetic_leapfrog)
        allocate(particle_kinetic_leapfrog::particles_tmp(blocksize))
      type is (particle_gc)
        allocate(particle_gc::particles_tmp(blocksize))
      type is (particle_gc_vpar)
        allocate(particle_gc_vpar::particles_tmp(blocksize))
      class default
        write(*,*) "ERROR: particle type not supported yet for initialize_particles_H_mu_psi"
        call exit(1)
    end select
    !allocate(particles_tmp(blocksize), mold=particles) ! this does not work in ifort 17
    allocate(found(blocksize))

#ifndef __NVCOMPILER
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
    !$omp   shared(particles_tmp, psimax, psimin, found, F0, cor, mass, charge, T_Maxwell, &
    !$omp          fields, psi_minmax_list, rans, R_axis, Z_axis, blocksize, &
    !$omp          my_include_vpar, central_density, init_uniform_space, Rbox, Zbox, space_pdf, use_rej) &
#endif
    !$omp   private(i, psi, theta, phi, i_elm, s, t, R, Z, R_s, R_t, Z_s, Z_t, &
    !$omp           u, particle_kinetic_tmp, v2, v_par,                        &
#ifdef fullmhd
    !$omp           A3, AR, AZ, A3_R, A3_Z, AR_Z, AR_p, AZ_R, AZ_P, Fprof,          &
#endif
    !$omp           P, P_s, P_t, P_phi, inv_st_jac, psi_R, psi_Z, B, H, muB, chi, ran, particle, temp, ifail, DUMMY_R, DUMMY_Z)
#endif
    do i=1,blocksize

      ran(:) = rans(:,i)

      if (init_uniform_space) then
        call transform_uniform_cylindrical([ran(3),ran(4),ran(5)], Rbox, Zbox, [0.d0,TWOPI], R, Z, phi)
        call find_RZ(fields%node_list, fields%element_list,R,Z,DUMMY_R,DUMMY_Z,i_elm,s,t,ifail)

        !> Reject this position if it fails the spatial pdf
        if (use_rej .and. i_elm .gt. 0) then
          if (eval_rej_f(fields%node_list, fields%element_list, i_elm, s, t, phi, R, Z, space_pdf) .lt. ran(7)) i_elm = 0
        end if

      else

        if (present(Psi_transform)) then
          psi = Psi_transform(ran(3))
        else
          psi = (psimax-psimin)*ran(3)+psimin
        end if

        ! Try to find this position
        if (present(Theta_transform)) then
          theta = Theta_transform(ran(4))
        else
          theta = TWOPI*ran(4)
        end if

        phi = TWOPI*ran(5)
        ! 1. Find R, Z corresponding to psi, theta
        call find_theta_psi(fields%node_list, fields%element_list,psi_minmax_list,theta,psi,phi,R_axis,Z_axis,i_elm,s,t,R,Z)
      end if
      ! By this point i_elm, s, t, R, Z, phi must be set

      if (i_elm .gt. 0) then
        found(i) = .true.
        ! If we are here a suitable position has been found
        particle%weight = 1.d0 ! This is needed because the initializing values for
        ! the particles are not being used always!
        particle%i_elm = i_elm
        particle%st    = [s,t]
        particle%x     = [R,Z,phi]

        ! 1. Get B at this position
#ifdef fullmhd
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[var_A3,var_AR,var_AZ],3,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
        call fields%calc_F_profile(i_elm,s,t,phi,Fprof)
        inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
        A3=P(1)
        AR=P(2)
        AZ=P(3)
        !Derivatives of A3
        A3_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
        A3_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac

        !Derivatives of AR
        AR_Z    = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
        AR_p    = P_phi(2)

        !Derivatives of AZ
        AZ_R    = (  P_s(3) * Z_t - P_t(3) * Z_s ) * inv_st_jac
        AZ_p    = P_phi(3)
        B=[(A3_Z-AZ_p)/R, (AR_p-A3_R)/R, AZ_R-AR_Z + Fprof/R]
        ! 2. Calculate E and mu, save in particle
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[var_T],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)

#else
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[1],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
        inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
        psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
        psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
        ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
        B        = [+psi_Z, -psi_R, F0] / R

        ! 2. Calculate E and mu, save in particle
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[6],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)


#endif
        if (present(T_Maxwell)) then
          temp = T_Maxwell
        else
          temp = P(1)/(2.d0*MU_ZERO*central_density*1.d20*EL_CHG) ! [eV]
#ifdef WITH_TiTe
          temp = temp*2d0 ! P(1) contains the ion temperature in this model, reverse previous correction
#endif
        endif

        H = temp*0.5d0*sample_chi_squared_3(ran(1))
        ! Solve now for u = 1-sqrt(1-x) (CDF of beta(1,1/2) distribution)
        ! Inverse gives: 1-(1-u**2) = x or x = 2u - u**2
        ! mu*B = v_par and mu*B/H = v_per**2/(v_per**2+v_par**2) ~ beta(1,1/2); where H is the total kinetic energy.
        ! ran(2) must be modified to randomize velocity properly, stored in u
        ! u and the sign for muB both need to have the same ran(2)
        u = 2*mod(ran(2),0.5d0)
        muB = sign(H*(2.d0*u-u**2), ran(2)-0.5d0)
        if (my_include_vpar) then
          call interp_PRZ(fields%node_list, fields%element_list,i_elm,[7],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
          ! Convert to sqrt of parallel energy
          P(1) = P(1)*sqrt(mass*ATOMIC_MASS_UNIT/EL_CHG)
          H = H + P(1)*(P(1) + 2.d0*sign(sqrt(H-muB),muB))
        end if

        particle%E  = H
        particle%mu = muB/norm2(B)

        ! 3. Calculate charge (if cor is present)
        if (present(cor)) then
          particle%q = int(q_coronal(fields%node_list, fields%element_list, s, t, phi, i_elm, cor, ran(8)),1)
        else
          if (present(charge)) then
            particle%q = int(charge,1)
          else
            particle%q = int(1,1) ! default value, warn here maybe?
          end if
        end if

        ! 4. Output to particles (dependent on type of particle)
        chi = TWOPI*ran(6)
        select type(p1 => particles_tmp(i))
          type is (particle_kinetic_leapfrog)

! the generic copy of particle_kinetic_leapfrog, i.e p = ..., seems broken, therefor using the non-generic copy
!          call copy_particle_kinetic_leapfrog( &
!                 kinetic_to_kinetic_leapfrog(gc_to_kinetic(fields%node_list, fields%element_list, particle, chi, B, mass), &
!                                             [0.d0, 0.d0, 0.d0], B, mass, dt=0.d0), &
!                                             p )
            particles_tmp(i) = gc_to_kinetic_leapfrog(particle, fields%node_list, fields%element_list, chi, [0.d0,0.d0,0.d0], B, mass, dt=0.d0)

            ! if the kinetic position is not in the grid particles(i)%i_elm the particle is lost
            if (particles_tmp(i)%i_elm .le. 0) found(i) = .false.
          type is (particle_gc)
            particles_tmp(i) = particle
          type is (particle_gc_vpar)
            call convert_gc_to_gc_vpar(particle, norm2(B), mass, p1)
        end select
      else
        found(i) = .false.
      end if
    end do
#ifndef __NVCOMPILER
    !$omp end parallel do
#endif
    ! How many particles have we found?
    n_found = count(found)
    write(*,*) my_id, "tried to find ", to_find, " found: ", n_found

    i=1
    do j=1,size(particles_tmp)
      if (found(j)) then
        if (i .gt. to_find) exit
        select type(p1 => particles(particles_done_local+i))
          type is (particle_kinetic_leapfrog)
            select type (p2 => particles_tmp(j))
              type is (particle_kinetic_leapfrog)
                p1 = p2
            end select
          type is (particle_gc)
            select type (p2 => particles_tmp(j))
              type is (particle_gc)
                p1 = p2
            end select
          type is (particle_gc_vpar)
            select type (p2 => particles_tmp(j))
            type is (particle_gc_vpar)
              p1 = p2
            end select
        end select

        i = i+1
      end if
    end do
    particles_done_local = particles_done_local + i-1

    deallocate(rans, found, particles_tmp)

    ! check if everyone is done
    call MPI_AllReduce(particles_to_do_local .eq. particles_done_local, all_done, &
      1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  end do
end subroutine initialise_particles_H_mu_psi

!> Subroutine for initialising particles in  mu, (psi, theta|R, Z), phi, gamma (gyrophase) space.
!! Optionally either only at a uniformly distributed n_phi_planes_in set of planes or with
!! multiple marker particles per guiding centre particle.
subroutine initialise_particles_H_mu_psi_phiplanes(particles, fields, rng_base, mass, T_maxwell, &
  Theta_transform, Psi_transform, alpha, E_max, include_vpar, uniform_space, &
  space_pdf, cor, charge, n_phi_planes_in, &
  n_gyro_orbit_in, rng_n_streams_round_off_in)
  use mod_rng
  use mod_fields
  use mod_random_seed
  use mod_sampling
  use constants
  use phys_module, only: F0, central_density
  use mod_coronal
  use mod_boris, only: gc_to_kinetic_leapfrog, kinetic_to_kinetic_leapfrog, gc_to_kinetic
  use mod_gc_variational, only: convert_gc_to_gc_vpar
  use mpi
  use mod_interp
  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  class(fields_base),    intent(in)                 :: fields
  class(type_rng),       intent(in)                 :: rng_base !< What type of random number generator to use (will be reseeded here)
  real*8,                intent(in)                 :: mass
  real*8,                external,   optional       :: Theta_transform !< Function to transform 0-1 to the theta-domain
  real*8,                external,   optional       :: Psi_transform !< Function to transform 0-1 to the Psi-domain
  !< if om               itted, determine auttically from node_list
  real*8,                intent(in), optional       :: alpha !< Make more fast (>0) or slow (<0) particles and weigh them appropriately
  real*8,                intent(in), optional       :: E_max !< If alpha=1 we select particles from a block-distribution, up to E_max
  logical,               intent(in), optional       :: include_vpar !< Initialize particles with local parallel velocity
  logical,               intent(in), optional       :: uniform_space !< Do not
  !< use {psi,theta}_transform if present bute rejection sampling in RZ
  type(spatial_pdf),     intent(in), optional       :: space_pdf !< spatial pdf to reject against when uniform_space is set. if absent the RZ sampling is uniform
  type(coronal),         intent(in), optional       :: cor !< Coronal equilibrium datatype for this particle. If unset, do not alter q
  integer,               intent(in), optional       :: charge !< Use this if cor is not present
  real*8,                intent(in), optional       :: T_Maxwell !< constant Maxwellian temperature [eV]
  integer,               intent(in), optional       :: n_phi_planes_in !Amount of phi planes
  integer,               intent(in), optional       :: n_gyro_orbit_in !Amount of phi planes
  logical,               intent(in), optional       :: rng_n_streams_round_off_in !< round-off the rng n_streams at 2**ceil

  ! Internal variables
  type(particle_gc) :: particle
  type(particle_kinetic) :: particle_kinetic_tmp
  class(type_rng), allocatable :: rng
  real*8  :: ran(8)
  real*8  :: H, muB, chi
  real*8  :: psi, psimin, psimax, theta, phi
  real*8  :: R, Z, inv_st_jac, psi_r, psi_z, B(3)
#ifdef fullmhd
  real*8    :: A3, AR, AZ, A3_R, A3_Z, AR_Z, AR_p, AZ_R, AZ_P, Fprof
  real*8, dimension(3)                :: P, P_s, P_t, P_phi
#else 
  real*8, dimension(1)                :: P, P_s, P_t, P_phi
#endif
  real*8  :: R_s, R_t, Z_s, Z_t
  real*8  :: s, t, u_init_max, temp, u, v2, v_par
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
  integer :: i_elm, i, j, ifail, my_id, n_mpi, ierr
  real*8, dimension(fields%element_list%n_elements,2)    :: psi_minmax_list
  real*8, allocatable, dimension(:,:)             :: rans
  class(particle_base), dimension(:), allocatable :: particles_tmp
  logical, dimension(:), allocatable              :: found
  real*8  :: t0, t1
  real*8  :: Rbox(2), Zbox(2), DUMMY_R, DUMMY_Z
  integer :: blocksize, prev_blocksize, particles_to_do_local, particles_done_local, blocksize_tmp
  integer :: to_find, n_tries_now, n_found
  logical :: all_done, init_uniform_space, my_include_vpar
  logical :: init_phiplanes, init_gyro_orbit, rng_n_streams_round_off, use_rej
  real*8  :: my_alpha
  integer :: n_phi_planes,i_phi_planes, n_gyro_orbit, i_gyro_temp, i_gyro_orbit

  rng_n_streams_round_off = .false.
  if(present(rng_n_streams_round_off_in)) rng_n_streams_round_off = rng_n_streams_round_off_in

  init_uniform_space = .false.
  if (present(uniform_space) .and. uniform_space) then
    init_uniform_space = .true.
    call domain_bounding_box(fields%node_list, fields%element_list, Rbox(1), Rbox(2), Zbox(1), Zbox(2))
  end if

  my_include_vpar = .false.
  if (present(include_vpar)) my_include_vpar = include_vpar
  ! Check if we are in the 3,4,5 series of models
  if (my_include_vpar .and. .not. with_Vpar) then
    write(*,*) "WARNING: This model, ", jorek_model, "does not support parallel flows, disabling"
    my_include_vpar = .false.
  end if

  !Phi planes
  init_phiplanes= .false.
  if (present(n_phi_planes_in)) then
    n_phi_planes=n_phi_planes_in
    init_phiplanes = .true.
    write(*,"(A,I4,A)") "WARNING: Initializing on ",n_phi_planes," planes around the torus"
  else
    n_phi_planes=1 !In this case, the loops aren't there
  endif
  !Gyro orbits
  init_gyro_orbit=.false.
  if (present(n_gyro_orbit_in)) then
    n_gyro_orbit=n_gyro_orbit_in
    init_gyro_orbit = .true.
    write(*,"(A,I4,A)") "WARNING: Initializing ",n_gyro_orbit," particles per guiding centre"
  else
    n_gyro_orbit=1 !In this case, the loops aren't there
  endif

 use_rej = .false.
 if (present(space_pdf)) then
  call check_spatial_pdf(space_pdf, "initialise_particles_H_mu_psi_phiplanes")
  use_rej = associated(space_pdf%f)
 endif

  if (present(alpha)) then
    write(*,*) "alpha not implemented yet"
    call exit(1)
    if (alpha .eq. 1.d0 .and. .not. present(E_max)) then
      write(*,*) "E_max is required if alpha=1"
      call exit(1)
    end if
    my_alpha = alpha
  else
    my_alpha = 0.d0
  end if

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ierr)

  if (.not. init_uniform_space) then
    psimin= 1d10
    psimax=-1d10
    ! Preparatory work: determine psi_min,max
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
#endif
    !$omp shared(fields,  psi_minmax_list) &
    !$omp private(i_elm) reduction(min:psimin) &
    !$omp reduction(max:psimax)
    do i_elm=1, fields%element_list%n_elements
      call psi_minmax(fields%node_list, fields%element_list, i_elm, psi_minmax_list(i_elm,1), psi_minmax_list(i_elm,2))
      psimin = min(psi_minmax_list(i_elm,1),psimin)
      psimax = max(psi_minmax_list(i_elm,2),psimax)
    end do
    !$omp end parallel do
  end if

  ! Preparatory work: setup RNG
  allocate(rng,source=rng_base)
  call rng%initialize(8, random_seed(), 1, 1, ifail, round_off_n_streams_in=rng_n_streams_round_off)

  ! Preparatory work: get R_axis, Z_axis
  if (.not. init_uniform_space) then
    call find_axis(my_id, fields%node_list, fields%element_list, psi_axis, R_axis, Z_axis, i_elm, s_axis, t_axis, ifail)
  end if

  ! Preset i_elm to 0 so that by default particles are lost
  particles(:)%i_elm = 0

  ! Assume equal number of particles everywhere
  particles_to_do_local  = size(particles)!/n_phi_planes !this is the local version
  particles_done_local   = 0
  prev_blocksize = 0 ! initial block size
  all_done = .false.
  !particles are stored locally

  do while (.not. all_done)
    to_find = (particles_to_do_local-particles_done_local)
    n_tries_now = to_find
    ! Add a little bit extra to cut off tail of 1/x^n
    if (n_tries_now .lt. 64 .and. n_tries_now .gt. 0) n_tries_now = 64

    ! Calculate starting-point for this block
    ! random numbers are generated in blocks per mpi process
    ! blocks are grouped in a superblock of block*n_mpi size
    ! blocks must be equal size
    ! we communicate a blocksize here, and keep a running index of the previous starting points
    ! make this the biggest of n_tries_now
    call MPI_AllReduce(n_tries_now, blocksize, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    call rng%jump_ahead((n_mpi-my_id-1)*prev_blocksize + my_id*blocksize)

    prev_blocksize = blocksize

    ! Generate the random numbers
    call cpu_time(t0)
    allocate(rans(8,blocksize))
    do i=1,blocksize
      call rng%next(rans(:,i))
    end do
    call cpu_time(t1)

    ! Allocate some temporary storage
    select type (particles)
      type is (particle_kinetic)
        write(*,*) "ERROR: particle_kinetic not supported yet for initialize_particles_H_mu_psi"
      type is (particle_kinetic_leapfrog)
        allocate(particle_kinetic_leapfrog::particles_tmp(blocksize))
      type is (particle_gc)
        allocate(particle_gc::particles_tmp(blocksize))
      class default
        write(*,*) "ERROR: particle type not supported yet for initialize_particles_H_mu_psi"
        call exit(1)
    end select
    !allocate(particles_tmp(blocksize), mold=particles) ! this does not work in ifort 17
    allocate(found(blocksize))
    !If we are initializing with gyro orbits, use reduced blocksize.
    if(init_gyro_orbit) then
      blocksize_tmp=floor(float(blocksize)/float(n_gyro_orbit))
    else
      blocksize_tmp=blocksize
    endif
    
    !If initializatin on a set of phi planes, make blocksize a multiple of the amount of phi planes
    if(init_phiplanes) then 
      blocksize_tmp=blocksize_tmp-modulo(blocksize_tmp,n_phi_planes)
      !Subtract as to keep within the bound of the blocksize array
    endif 
#ifndef __NVCOMPILER
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
    !$omp   shared(particles_tmp, psimax, psimin, found, F0, cor, mass, charge, T_Maxwell, &
    !$omp          fields, psi_minmax_list, rans, R_axis, Z_axis, blocksize, init_phiplanes,init_gyro_orbit, n_gyro_orbit,blocksize_tmp,&
    !$omp          my_include_vpar, central_density, init_uniform_space, Rbox, Zbox, space_pdf, use_rej, n_phi_planes,my_id) &
#endif
    !$omp   private(i, psi, theta, phi, i_elm, s, t, R, Z, R_s, R_t, Z_s, Z_t, &
    !$omp           u, particle_kinetic_tmp, v2, v_par,                        &
#ifdef fullmhd
    !$omp          A3, AR, AZ, A3_R, A3_Z, AR_Z, AR_p, AZ_R, AZ_P, Fprof, &
#endif
    !$omp           P, P_s, P_t, P_phi, inv_st_jac, psi_R, psi_Z, B, H, muB, chi, ran, particle, temp, ifail, DUMMY_R, DUMMY_Z,i_phi_planes,&
    !$omp           i_gyro_orbit, i_gyro_temp)
#endif
    do i=1,blocksize_tmp

      ran(:) = rans(:,i)

      if (init_uniform_space) then
        if( init_phiplanes) then
          !Plane initialization, all particles sampled on the 0 plane.
          call transform_uniform_cylindrical([ran(3),ran(4),ran(5)], Rbox, Zbox, [0.d0,0.d0], R, Z, phi)
        else
          call transform_uniform_cylindrical([ran(3),ran(4),ran(5)], Rbox, Zbox, [0.d0,TWOPI], R, Z, phi)
        endif


        call find_RZ(fields%node_list, fields%element_list,R,Z,DUMMY_R,DUMMY_Z,i_elm,s,t,ifail)

        !> reject this position if it fails the spatial pdf
        if (use_rej .and. i_elm .gt. 0) then
          if (eval_rej_f(fields%node_list, fields%element_list, i_elm, s, t, phi, R, Z, space_pdf) .lt. ran(7)) i_elm = 0
        end if

      else

        if (present(Psi_transform)) then
          psi = Psi_transform(ran(3))
        else
          psi = (psimax-psimin)*ran(3)+psimin
        end if

        ! Try to find this position
        if (present(Theta_transform)) then
          theta = Theta_transform(ran(4))
        else
          theta = TWOPI*ran(4)
        end if

        phi =TWOPI*ran(5)
        ! 1. Find R, Z corresponding to psi, theta
        call find_theta_psi(fields%node_list, fields%element_list,psi_minmax_list,theta,psi,phi,R_axis,Z_axis,i_elm,s,t,R,Z)
      end if
      ! By this point i_elm, s, t, R, Z, phi must be set

      if (i_elm .gt. 0) then
        !Can only set found later
        ! If we are here a suitable position has been found
        particle%weight = 1. ! This is needed because the initializing values for
        ! the particles are not being used always!
        particle%i_elm = i_elm
        particle%st    = [s,t]
        particle%x     = [R,Z,phi]

        ! 1. Get B at this position

        !Changed for Full MHD. Note that in equilibrium, these are basically equivalent (no AR & AZ). If initializing on non-equilibrium
        !fluid states, these corrections are needed.
#ifdef fullmhd
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[var_A3,var_AR,var_AZ],3,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
        call fields%calc_F_profile(i_elm,s,t,phi,Fprof)
        inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
        A3=P(1)
        AR=P(2)
        AZ=P(3)
        !Derivatives of A3
        A3_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
        A3_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac

        !Derivatives of AR
        AR_Z    = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
        AR_p    = P_phi(2)

        !Derivatives of AZ
        AZ_R    = (  P_s(3) * Z_t - P_t(3) * Z_s ) * inv_st_jac
        AZ_p    = P_phi(3)
        B=[(A3_Z-AZ_p)/R, (AR_p-A3_R)/R, AZ_R-AR_Z + Fprof/R]
        ! 2. Calculate E and mu, save in particle
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[var_T],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
#else
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[1],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
        inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
        psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
        psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
        ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
        B        = [+psi_Z, -psi_R, F0] / R

        ! 2. Calculate E and mu, save in particle
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[6],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)

#endif
        if (present(T_Maxwell)) then
          temp = T_Maxwell
        else
          temp = P(1)/(2.d0*MU_ZERO*central_density*1.d20*EL_CHG) ! [eV]
#ifdef WITH_TiTe
          temp = temp*2d0 ! P(1) contains the ion temperature in this model, reverse previous correction
#endif
        endif

        H = temp*0.5d0*sample_chi_squared_3(ran(1))
        ! Solve now for u = 1-sqrt(1-x) (CDF of beta(1,1/2) distribution)
        ! Inverse gives: 1-(1-u**2) = x or x = 2u - u**2
        ! mu*B = v_par and mu*B/H = v_per**2/(v_per**2+v_par**2) ~ beta(1,1/2); where H is the total kinetic energy.
        ! ran(2) must be modified to randomize velocity properly, stored in u
        ! u and the sign for muB both need to have the same ran(2)
        u = 2*mod(ran(2),0.5d0)
        muB = sign(H*(2.d0*u-u**2), ran(2)-0.5d0)
        if (my_include_vpar) then
          call interp_PRZ(fields%node_list, fields%element_list,i_elm,[7],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
          ! Convert to sqrt of parallel energy
          P(1) = P(1)*sqrt(mass*ATOMIC_MASS_UNIT/EL_CHG)
          H = H + P(1)*(P(1) + 2.d0*sign(sqrt(H-muB),muB))
        end if

        particle%E  = H
        particle%mu = muB/norm2(B)

        ! 3. Calculate charge (if cor is present)
        if (present(cor)) then
          particle%q = int(q_coronal(fields%node_list, fields%element_list, s, t, phi, i_elm, cor, ran(8)),1)
        else
          if (present(charge)) then
            particle%q = int(charge,1)
          else
            particle%q = int(1,1) ! default value, warn here maybe?
          end if
        end if

        ! 4. Output to particles (dependent on type of particle)
        chi = TWOPI*ran(6)
        do i_gyro_orbit=1,n_gyro_orbit
          i_gyro_temp=(i-1)*n_gyro_orbit+i_gyro_orbit

          if (i*n_gyro_orbit .gt. blocksize) exit
          select type(p1 => particles_tmp(i_gyro_temp))
            type is (particle_kinetic_leapfrog)

              ! the generic copy of particle_kinetic_leapfrog, i.e p = ..., seems broken, therefore using the non-generic copy

!              call copy_particle_kinetic_leapfrog( &
!                kinetic_to_kinetic_leapfrog(gc_to_kinetic(fields%node_list, fields%element_list, particle, TWOPI*REAL(i_gyro_orbit)/REAL(n_gyro_orbit)+chi, B, mass), &
!                [0.d0, 0.d0, 0.d0], B, mass, dt=0.d0), &
!                p )

              particles_tmp(i_gyro_temp) = gc_to_kinetic_leapfrog(particle, fields%node_list, fields%element_list, TWOPI*REAL(i_gyro_orbit,8)/REAL(n_gyro_orbit,8)+chi, [0.d0,0.d0,0.d0], B, mass, dt=0.d0)
  
              ! if the kinetic position is not in the grid particles(i)%i_elm the particle is lost
              found(i_gyro_temp) = .true.
              if (p1%i_elm .le. 0) found(i_gyro_temp) = .false.

            type is (particle_gc)
              p1 = particle
            type is (particle_gc_vpar)
              call convert_gc_to_gc_vpar(particle, norm2(B), mass, p1)
          end select !particle type
        end do !n_gyro_orbit

      else
        do i_gyro_orbit=1,n_gyro_orbit
          found((i-1)*n_gyro_orbit+i_gyro_orbit) = .false. !if gyro_orbit initialized, need to set a range as not found
        enddo
      end if
    end do
#ifndef __NVCOMPILER
    !$omp end parallel do
#endif

    ! How many particles have we found?
    n_found = count(found)
    write(*,*) my_id, "tried to find ", to_find, " found: ", n_found

    i=1
    do j=1,size(particles_tmp)
      if (found(j)) then
        do i_phi_planes=1,n_phi_planes
          if (i .gt. to_find) exit
          select type(p1 => particles(particles_done_local+i))
            type is (particle_kinetic_leapfrog)
              select type (p2 => particles_tmp(j))
                type is (particle_kinetic_leapfrog)
                  if (init_phiplanes) then
                    p2%x(3)=p2%x(3)+TWOPI/REAL(n_phi_planes) !+= because particles are tilted wrt poloidal plane due to B != B_phi
                    !p2 is modified by this staement, keeps pointing to the same particle -> add 2\pi /n_phi_planes every time.
                  endif
                  p1 = p2
              end select
            type is (particle_gc)
              select type (p2 => particles_tmp(j))
                type is (particle_gc)
                  if (init_phiplanes) then
                    p2%x(3)=p2%x(3)+TWOPI/REAL(n_phi_planes)
                  endif
                  p1 = p2
              end select
          end select
          i = i+1
        enddo !phi planes
      end if !found particle
    end do !particles in particle_tmp

    particles_done_local = particles_done_local + i-1

    deallocate(rans, found, particles_tmp)

    ! check if everyone is done
    call MPI_AllReduce(particles_to_do_local .eq. particles_done_local, all_done, &
      1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
  end do

end subroutine initialise_particles_H_mu_psi_phiplanes

!> Calculate the particle weights according to the canonical maxwellian distribution function
!> (no electric fields):
!> \[ F_{MC}(P_\phi,H,\mu) = \frac{n(\bar\psi)}{\left[2\pi \bar T(\bar\psi)/m\right]^3/2}
!>                          exp\left{-\frac{H}{\bar T(\bar\psi)}\right} \]
!> where \(\bar\psi = P_\phi/q\), \(P_\phi = q\psi + m R v_\phi\),
!> \( H = m/2 v_\parallel^2 + \mu B \) and \(\mu = \frac{m v_\perp^2}{2B}\)
!>
!> The particle weight is set to the value of this distribution function at the specific point.
!> \(\bar T(\bar\psi)\) is approximated by \(T(\psi)\) if it is missing.
!> \(\bar n(\bar\psi)\) is approximated by 1 if it is missing.
subroutine set_particle_weights_canonical_maxwellian(particles, node_list, element_list, mass, n_psibar, T_psibar, alpha)
  use data_structure
  use constants
  use phys_module, only: central_density, central_mass
  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  real*8, intent(in)                                :: mass
  real*8, external, optional                        :: n_psibar
  real*8, external, optional                        :: T_psibar !< Provide the temperature in eV as a function of psibar
  real*8, optional                                  :: alpha !< Interpolation factor between sampling from a block (1) and the maxwellian (0)

  integer              :: i
  real*8               :: t_norm, psibar, H, n, T
  real*8, dimension(1) :: P, P_s, P_t, P_phi
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t
  real*8               :: my_alpha

  t_norm = sqrt(MU_ZERO * central_mass * ATOMIC_MASS_UNIT * central_density * 1.d20)

  if (present(alpha)) then
    my_alpha = alpha
  else
    my_alpha = 0.d0 !< Default
  end if

  ! default(shared) is very dangerous but needed due to gfortran failures... be careful adding variables
  ! and try compilation with default(none) if you change anything.
#ifdef __GFORTRAN__
  !$omp parallel do default(shared) &
#else
  !$omp parallel do default(none) &
#endif
  !$omp private(i, psibar, H, n, T, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t) &
  !$omp shared(particles, node_list, element_list, mass, central_density, my_alpha)
  do i=1,size(particles,1)

    if (particles(i)%i_elm .le. 0) cycle

    call interp_PRZ(node_list,element_list,particles(i)%i_elm,[1],1,        &
      particles(i)%st(1),particles(i)%st(2),particles(i)%x(3), &
      P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)

    select type (pa => particles(i))
      type is (particle_kinetic_leapfrog)
        psibar = real(pa%q,8) * P(1) * EL_CHG + mass * ATOMIC_MASS_UNIT * R * pa%v(3)
        H = mass*ATOMIC_MASS_UNIT*0.5d0 * norm2(pa%v) / EL_CHG ! [eV]
      class default
        write(*,*) "ERROR: add code for your type here"
        cycle
    end select

    if (present(n_psibar)) then
      n = n_psibar(psibar) ! [m^-3]
    else
      n = 1 ! units irrelevant if normalized later to a total number of particles
    end if

    if (present(T_psibar)) then
      T = T_psibar(psibar) ! [eV]
    else
      ! Calculate the local temperature and use this instead
      call       interp_PRZ(node_list,element_list,particles(i)%i_elm,[6],1, &
        particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
      ! P(1)/(kb mu_zero n_zero) is in [K], multiply by kb/el_chg to go to eV
      T = P(1)/(2.d0*MU_ZERO*central_density*1.d20*EL_CHG) ! [eV] factor 2 is due to
      ! definition of P(1) as ion + electron temperature
#ifdef WITH_TiTe
      T = T*2d0 ! P(1) contains the ion temperature in this model, reverse previous correction
#endif
      ! Workaround for low-temperature regions
      ! Because we give weights based on the energy and temperature areas with lower temperature are getting
      ! too high weights. Work around this by defining a minimum temperature to stop the outer regions from
      ! dominating the projection.
      T = max(1d1,T)
    end if
    particles(i)%weight = real(n * exp(-H*my_alpha/T),8) ! if zero, all particles have equal weight n
  end do
  !$omp end parallel do
end subroutine set_particle_weights_canonical_maxwellian

!> Normalize particles with the result of the projection of the first group
subroutine normalize_with_projection(proj, particles, i_group)
  use mod_project_particles
  use mod_interp, only: interp_0
  type(projection), intent(in) :: proj
  class(particle_base), dimension(:), intent(inout) :: particles
  integer, intent(in), optional :: i_group

  integer :: group = 1
  integer :: i
  real*8, dimension(1) :: P
  if (present(i_group)) group = i_group

  do i=1,size(particles,1)
    if (particles(i)%i_elm .gt. 0) then
      call interp_0(proj%node_list,proj%element_list,particles(i)%i_elm,[group],1, &
        particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),P)
      P = max(P, 1d-2) ! guard against divide-by-zero, maximum adjustment ratio is then 10^2
      particles(i)%weight = real(particles(i)%weight/P(1),8)
    end if
  end do
end subroutine normalize_with_projection

!> Normalize particles with the result of the projection of i_group
subroutine normalize_with_projection_at_gc(proj, particles, fields, time, mass, i_group)
  use mod_fields
  use mod_project_particles
  use mod_interp, only: interp_0
  use mod_boris,  only: kinetic_leapfrog_to_gc
  type(projection), intent(in)                      :: proj
  class(particle_base), dimension(:), intent(inout) :: particles
  class(fields_base), intent(in)                    :: fields

  real*8, intent(in)            :: time
  real*8, intent(in)            :: mass
  integer, intent(in), optional :: i_group

  integer              :: group = 1
  real*8               :: my_dt = 0.d0
  integer              :: i
  real*8, dimension(1) :: P
  real*8               :: B(3), E(3), psi, U
  type(particle_gc)    :: p_gc

  if (present(i_group)) group = i_group

  do i=1,size(particles,1)

    if (particles(i)%i_elm .gt. 0) then

      call fields%calc_EBpsiU(time, particles(i)%i_elm, particles(i)%st, particles(i)%x(3), E, B, psi, U)

      select type (p => particles(i))
        type is (particle_kinetic_leapfrog)
          p_gc = kinetic_leapfrog_to_gc(fields%node_list, fields%element_list, p, [0.d0,0.d0,0.d0], B, mass, dt=0.d0)
        type is (particle_gc)
          p_gc = p
      end select

      if (p_gc%i_elm .gt. 0) then
        call interp_0(proj%node_list,proj%element_list,p_gc%i_elm,[group],1, p_gc%st(1),p_gc%st(2),p_gc%x(3),P)
        P = max(P, 1d-2) ! guard against divide-by-zero, maximum adjustment ratio is then 10^2
        particles(i)%weight = real(particles(i)%weight/P(1),8)
      else
        call interp_0(proj%node_list,proj%element_list,particles(i)%i_elm,[group],1, &
          particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),P)
        P = max(P, 1d-2) ! guard against divide-by-zero, maximum adjustment ratio is then 10^2
        particles(i)%weight = real(particles(i)%weight/P(1),8)
        ! otherwise weigh at normal position to not screw up the weighting
      end if
    end if
  end do
end subroutine normalize_with_projection_at_gc


!> Weigh particles with interpolated values.
!> Note that this multiplies existing weights by the value of f, instead of
!> overwriting weights. This allows using it with nonuniform weights.
subroutine weigh_with_interp_f(node_list, element_list, particles, vars, f)
  use mod_interp
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  class(particle_base), dimension(:), intent(inout) :: particles
  integer, dimension(:), intent(in)                 :: vars  !< Which variables to interpolate
  !< Special numbers: -1: R, -2: Phi, -3: Z. Must be in ascending order!
  interface
    !> Transformation function on interpolated variables
    pure function f(P)
      real*8, intent(in), dimension(:) :: P
      real*4 :: f
    end function f
  end interface

  integer :: i, k, n_mhd, n_geom
  real*8, dimension(size(vars,1)) :: P

  ! Get the number of mhd variables to use
  n_mhd = count(vars .gt. 0)
  n_geom = size(vars, 1) - n_mhd

  do i=1,size(particles,1)

    if (particles(i)%i_elm .gt. 0) then

      call interp_0(node_list,element_list,particles(i)%i_elm,                             &
        vars(n_geom:n_geom+n_mhd),n_mhd,particles(i)%st(1),particles(i)%st(2), &
        particles(i)%x(3), P(n_geom:n_geom+n_mhd))

      do k=1,n_geom
        select case (vars(k))
          case (0);  P(k) = 1.d0
          case (-1); P(k) = particles(i)%x(1)
          case (-2); P(k) = particles(i)%x(2)
          case (-3); P(k) = particles(i)%x(3)
        end select
      end do
      particles(i)%weight = particles(i)%weight * f(P)

    end if

  end do
end subroutine weigh_with_interp_f


!> Calculate the size of a box around the domain (in RZ)
subroutine domain_bounding_box(node_list, element_list, Rmin, Rmax, Zmin, Zmax)
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  real*8, intent(out)                               :: Rmin, Rmax, Zmin, Zmax
  real*8 :: el_Rmin, el_Rmax, el_Zmin, el_Zmax
  integer :: i
  ! Initial setting
  call RZ_minmax(node_list, element_list, 1, Rmin, Rmax, Zmin, Zmax)
  do i=2,element_list%n_elements
    call RZ_minmax(node_list, element_list, i, el_Rmin, el_Rmax, el_Zmin, el_Zmax)
    Rmin = min(Rmin, el_Rmin)
    Rmax = max(Rmax, el_Rmax)
    Zmin = min(Zmin, el_Zmin)
    Zmax = max(Zmax, el_Zmax)
  enddo
end subroutine domain_bounding_box


!> Adjust weights on all particles to have the correct number of atoms in total
subroutine adjust_particle_weights(particles, num_atoms_total)
  use mpi
  class(particle_base), intent(inout), dimension(:) :: particles
  real*8, intent(in)                                :: num_atoms_total !< What the sum of the weights should be
  real*8 :: local_weights, sum_weights, local_weights_active, sum_weights_active
  integer :: ifail, i

!  local_weights = sum(particles(:)%weight)
!  local_weights_active = sum(particles(:)%weight, mask=particles(:)%i_elm .gt. 0)

  local_weights        = 0d0
  local_weights_active = 0d0
  do i=1,size(particles)
    local_weights = local_weights + particles(i)%weight
    if (particles(i)%i_elm .gt. 0) local_weights_active = local_weights_active + particles(i)%weight
  enddo

  call MPI_AllReduce(local_weights,       sum_weights,       1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ifail)
  call MPI_AllReduce(local_weights_active,sum_weights_active,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ifail)

! Divide all weights by the sum of weights and multiply by the requested number of atoms
!  particles(:)%weight = real(real(particles(:)%weight,8) / sum_weights * num_atoms_total,4)

!  particles(:)%weight = particles(:)%weight / sum_weights_active * num_atoms_total

  do i=1, size(particles)
    particles(i)%weight = particles(i)%weight / sum_weights_active * num_atoms_total
  enddo

end subroutine adjust_particle_weights

function q_coronal(node_list, element_list, s, t, phi, i_elm, cor, u)
  use data_structure
  use phys_module, only: central_density
  use mod_coronal
  use mod_sampling, only: sample_discrete
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  real*8, intent(in)                                :: s, t, phi
  integer, intent(in)                               :: i_elm
  type(coronal), intent(in)                         :: cor
  integer                                           :: q_coronal
  real*8, intent(in)                                :: u !< uniform random number on [0,1]

  real*8, dimension(2) :: P, P_s, P_t, P_phi
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t
  real*8 :: local_Te, local_Ne, DUMMY_REAL
  real*8 :: q(0:cor%n_Z)

#ifdef WITH_TiTe
  call interp_PRZ(node_list,element_list,i_elm,[5,8],2,s,t,phi,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t) ! electron temperature
#else
  call interp_PRZ(node_list,element_list,i_elm,[5,6],2,s,t,phi,P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t) ! electron temperature + ion temperature (assumed equal)
#endif

  local_Ne = P(1) * 1d20                           ! plasma density [1/m^3]
  local_Te = P(2)/(2.d0*MU_ZERO*central_density*1.d20)/K_BOLTZ
#ifdef WITH_TiTe
  local_Te = local_Te*2.d0 ! P(1) contains the electron temperature, reverse previous correction
#endif

  if (local_Ne .le. 0.d0 .or. local_Te .le. 0.d0) then
    q_coronal = 0_1
  else
    call cor%interp(log10(local_Ne),log10(local_Te),q)
    q_coronal = int(sample_discrete(q, u),1)
  endif
end function

!> Set v of a particle for use with kinetic codes
subroutine set_velocity_from_T(particles, mass, node_list, element_list, cor, v_par)
  use constants
  use data_structure
  use phys_module, only: central_density, central_mass, F0
  use mod_sampling
  use mpi
  use mod_random_seed
  use mod_coordinate_transforms
  use mod_coronal
  implicit none

  class(particle_base), intent(inout), dimension(:) :: particles !< Particle to initialize
  real*8, intent(in)                                :: mass !< [u]
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  type(coronal), intent(in), optional               :: cor !< Coronal equilibrium datatype for this particle. If unset, do not alter q
  logical, intent(in), optional                     :: v_par !< Include the parallel velocity if present and true

  class(type_rng), allocatable :: my_rng
  integer :: i, ifail, seed, my_id, n_mpi
  real*8, dimension(4) :: P, P_s, P_t, P_phi
  real*8 :: v_out(3)
  real*8 :: R, R_s, R_t, Z, Z_s, Z_t, Psi, Psi_R, Psi_Z, B(3)
  real*8, parameter :: r_hat(3) = [1.d0, 0.d0, 0.d0]
  real*8 :: background_kbT, background_Kelvin, background_density, V_thermal
  real*8 :: DUMMY_REAL, t_norm
  real*8, allocatable :: Z_coronal(:)

  t_norm = sqrt(MU_ZERO * central_mass * ATOMIC_MASS_UNIT * central_density * 1.d20)

! Calculate a single random seed and communicate it over MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ifail)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ifail)
  if (my_id .eq. 0) seed = random_seed()
  call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ifail)

#ifndef WITH_Vpar
  if (present(v_par) .and. v_par) then
    write(*,*) "ERROR: initialization with v// not possible with this model"
    call MPI_ABORT(-1, MPI_COMM_WORLD, ifail)
  end if
#endif

  do i=1,size(particles)
    if (particles(i)%i_elm .le. 0) cycle
#ifdef WITH_TiTe
    call interp_PRZ(node_list,element_list,particles(i)%i_elm,[1,5,var_Te,7],4,particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),&
      P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
#else
    call interp_PRZ(node_list,element_list,particles(i)%i_elm,[1,5,var_T,7],4,particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),&
      P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
#endif

    background_density = P(2) * 1d20                           ! plasma density [1/m^3]
    ! Assume that the particles have the same temperature as the electrons
#ifdef WITH_TiTe
    background_kbT = P(3)/(MU_ZERO*central_density*1.d20)      ! P(1) contains the electron temperature
#else
    background_kbT = P(3)/(2.d0*MU_ZERO*central_density*1.d20) ! P(1) contains the total plasma temperature in J/kB = T = Te + Ti
#endif
    V_thermal = sqrt(background_kbT / (mass*ATOMIC_MASS_UNIT))      ! variance in each of the velocity dimensions [m/s]

    ! Only an implementation for particle_kinetic_leapfrog now
    select type (pa => particles(i))
      type is (particle_kinetic_leapfrog)

        ! v_out now contains parallel and perpendicular velocities and the gyrophase
        v_out(1:2) = boxmueller_transform(pa%v(1:2))*V_thermal ! 2 gaussian distributed random numbers
        v_out(3)   = sample_gaussian(pa%v(3))*V_thermal ! very slow, don't use in production

        if (present(cor)) then
          if (allocated(Z_coronal)) deallocate(Z_coronal)
          allocate(Z_coronal(0:cor%n_Z))
          background_kelvin  = background_kbT / K_BOLTZ              ! electron temperature [K]
          if (background_density .le. 0.d0 .or. background_kelvin .le. 0.d0) then
            Z_coronal = 0.d0
            Z_coronal(0) = 1.d0
          else
            call cor%interp(log10(background_density),log10(background_kelvin),Z_coronal)
          endif
        end if

        ! Calculate b^ (unit vector in direction of B)
        psi_R = (  P_s(4) * Z_t - P_t(4) * Z_s )/(R_s * Z_t - R_t * Z_s)
        psi_Z = (- P_s(4) * R_t + P_t(4) * R_s )/(R_s * Z_t - R_t * Z_s)
        B = [psi_Z, -psi_R, F0]/(R)

        ! Transform parallel and perpendicular velocities to R, Z, Phi
        ! To get the perpendicular vector, get a single vector perpendicular to b (b x r)
        ! and rotate it by another vector perpendicular to b.
        ! use the vector triple product to simplify.
        ! I'm not sure if this formula is the same in a right-handed coordinate system...
        ! this might change the direction of the rotation, but that is not important.
        if (present(v_par) .and. v_par) then
          pa%v = v_out(1) + (P(4)/t_norm) * B ! See normalisation of v_par
        else
          pa%v = v_out
        end if

        if (present(cor)) pa%q = int(maxloc(Z_coronal,1),1) ! take the most probable one here.
        ! should be better, with a random number and selection by probability
      class default
        write(*,*) "set_velocity_from_T not implemented for this particle type"
    end select
  end do
end subroutine set_velocity_from_T
end module initialisers_base
