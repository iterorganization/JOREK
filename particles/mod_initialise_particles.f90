!> Module for the initialization of particles in configuration space (6D)
!> by rejection sampling
module mod_initialise_particles
use mod_rng
use data_structure
use mod_particle_types
use constants
use mod_interp
implicit none
private
public initialise_particles, no_transform, adjust_particle_weights
public set_velocity_from_T, domain_bounding_box, initialise_particles_H_mu_psi
public set_particle_weights_canonical_maxwellian, normalize_with_projection
public weigh_with_interp_f
public normalize_with_projection_at_gc

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
  function rej_f(n, P, gradP)
    integer, intent(in) :: n
    real*8, dimension(n), intent(in) :: P
    real*8, dimension(3,n), intent(in) :: gradP
    real*4 :: rej_f
  end function rej_f
end interface
contains
!> Set positions for particles by rejection sampling from geometric and mhd
!> variables after collecting with transform, within Rbound, Zbound and Phibound
!> if present. See [[test_rejection_sampling]] for examples.
subroutine initialise_particles(particles, node_list, element_list, &
        rng, variables, transform, f, Rbound, Zbound, Phibound)
  use mpi
  use mod_sampling
  use mod_random_seed
  use mod_interp
  !$ use omp_lib
  implicit none

  class(particle_base), dimension(:), intent(inout) :: particles
  type(type_node_list), intent(in)                  :: node_list
  type(type_element_list), intent(in)               :: element_list
  class(type_rng), intent(in)                       :: rng !< What type of random number generator to use. Is re-seeded in the subroutine.
  integer, dimension(:), intent(in), optional       :: variables !< Which variables from JOREK to use. If absent, sample uniformly.
  real*8, external, optional                        :: transform !< Merge variables into a single criterium between 0 and 1 for rej.  sampling
  !< Special values: 0 = 1, -1 = R, -2 = Z, -3 = Phi. Must be in ascending order!
  real*8, intent(in), optional                      :: f !< Weighting factor: f=0 indicates uniform weights, f=1 indicates uniform distribution
  !< (particle weight proportional to transform(P) at that point.) If omitted take f=0.
  real*8, dimension(2), intent(in), optional        :: Rbound, Zbound, Phibound !< Between which coordinates to sample (RZPhi).
  !< if omitted, determine automatically from node_list

  ! Internal variables
  real*8  :: R, Z, phi, s, t, DUMMY_REAL
  real*8  :: Rbox(2), Zbox(2), Phibox(2)
  integer :: i, j, k, ifail
  real*8  :: ran(7)
  integer :: i_elm
  real*8  :: t0, t1, ostart, oend
  integer :: seq, n_streams, n_threads, i_thread
  integer :: n_geom, n_mhd
  integer :: my_id, n_cpu
  integer :: seed
  real*8, dimension(:), allocatable :: P
  class(type_rng), allocatable, dimension(:) :: rngs ! The RNGs for all the threads
  integer, dimension(:), allocatable :: i_to_find
  logical, dimension(:), allocatable :: not_found

  ostart = 0.d0
  oend   = 0.d0

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ifail)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ifail)
  if (present(variables)) then
    if (.not. present(transform)) then
      write(*,*) "ERROR: if variables are present in set_particle_position_rejection_sampling transform must also be present"
      call MPI_ABORT(MPI_COMM_WORLD, 10, ifail)
    end if
    ! Get the number of mhd variables to use
    allocate(P(size(variables,1)))
    n_mhd = count(variables .gt. 0)
    n_geom = size(variables, 1) - n_mhd
  else
    n_mhd = 0
    n_geom = 0
  end if

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
  n_streams = n_cpu*n_threads ! Works only for homogeneous environments!

  do i_thread=0,n_threads-1
    seq = my_id*n_threads + i_thread + 1
    call rngs(i_thread)%initialize(7, seed, n_streams, seq, ifail)
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
  !$omp   shared(particles, node_list, element_list, Rbox, Zbox, PhiBox, variables, &
  !$omp          rngs, n_threads, n_streams, seed, my_id, n_mhd, n_geom, i_to_find, not_found) &
  !$omp   private(j, i, R, Z, phi, i_elm, s, t, ifail, seq, ran, i_thread, P, DUMMY_REAL)
  i_thread = 0
  !$ i_thread=omp_get_thread_num()
  !$omp do schedule(static)
  do i=1,size(i_to_find,1)
    j = i_to_find(i)
    ! Generate a random position to put this particle
    call rngs(i_thread)%next(ran)
    call transform_uniform_cylindrical(ran(1:3), Rbox, Zbox, PhiBox, R, Z, phi)

    call find_RZ(node_list,element_list,R,Z,DUMMY_REAL,DUMMY_REAL,i_elm,s,t,ifail)
    if (ifail .eq. 0) then
      if (present(variables)) then
        ! Select the mhd variables requested
        if (n_mhd .ge. 1) then
          call interp_0(node_list,element_list,i_elm,variables(n_geom:n_geom+n_mhd),n_mhd,s,t,phi,P(n_geom:n_geom+n_mhd))
        end if
        do k=1,n_geom
          select case (variables(k))
          case (0);  P(k) = 1.d0
          case (-1); P(k) = R
          case (-2); P(k) = Z
          case (-3); P(k) = phi
          end select
        end do

        if (present(transform)) then
          if (ran(4) .lt. transform(p)) then
            particles(j)%x = [r, z, phi]
            particles(j)%i_elm = i_elm
            particles(j)%st = [s, t]
            select type (pa => particles(j))
            type is (particle_kinetic_leapfrog)
              pa%v = ran(5:7) ! save other components of this point for velocity init in a later routine
            end select
            not_found(i) = .false.
          end if
        end if
      else
        particles(j)%x = [r, z, phi]
        particles(j)%i_elm = i_elm
        particles(j)%st = [s, t]
        select type (pa => particles(j))
        type is (particle_kinetic_leapfrog)
          pa%v = ran(5:7) ! save other components of this point for velocity init in a later routine
        end select
        not_found(i) = .false.
      end if
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

!> Initialise particle positions in E, mu, (psi, theta|R, Z), phi, gamma (gyrophase) space.
!> Set Psi_transform to transform from [0,1] to your desired range
subroutine initialise_particles_H_mu_psi(particles, fields, rng_base, mass, T_maxwell, &
        Theta_transform, Psi_transform, alpha, E_max, include_vpar, uniform_space, &
        uniform_space_rej_f, uniform_space_rej_vars, cor, charge)
  use mod_rng
  use mod_fields
  use mod_random_seed
  use mod_sampling
  use constants
  use phys_module, only: F0, central_density
  use mod_coronal
  use mod_boris, only: gc_to_kinetic_leapfrog, kinetic_to_kinetic_leapfrog, gc_to_kinetic
  use mpi
  use mod_interp
  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  class(fields_base), intent(in)                    :: fields
  class(type_rng), intent(in)                       :: rng_base !< What type of random number generator to use (will be reseeded here)
  real*8, intent(in)                                :: mass
  real*8, external, optional                        :: Theta_transform !< Function to transform 0-1 to the theta-domain
  real*8, external, optional                        :: Psi_transform !< Function to transform 0-1 to the Psi-domain
  !< if omitted, determine automatically from node_list
  real*8, intent(in), optional                      :: alpha !< Make more fast (>0) or slow (<0) particles and weigh them appropriately
  real*8, intent(in), optional                      :: E_max !< If alpha=1 we select particles from a block-distribution, up to E_max
  logical, intent(in), optional                     :: include_vpar !< Initialize particles with local parallel velocity
  logical, intent(in), optional                     :: uniform_space !< Do not
  !< use {psi,theta}_transform if present but use rejection sampling in RZ
  procedure(rej_f), optional                        :: uniform_space_rej_f !< Merge variables into a single criterium between 0 and 1 for rej.  sampling
  !< Special values: 0 = 1, -1 = R, -2 = Z, -3 = Phi. Must be in ascending order!
  integer, dimension(:), intent(in), optional       :: uniform_space_rej_vars !< Variables to use for uniform_space_rej_f
  type(coronal), intent(in), optional               :: cor !< Coronal equilibrium datatype for this particle. If unset, do not alter q
  integer, intent(in), optional                     :: charge !< Use this if cor is not present
  real*8, intent(in), optional                      :: T_Maxwell !< constant Maxwellian temperature [eV]

  ! Internal variables
  type(particle_gc) :: particle
  class(type_rng), allocatable :: rng
  real*8  :: ran(8)
  real*8  :: H, muB, chi
  real*8  :: psi, psimin, psimax, theta, phi
  real*8  :: R, Z, inv_st_jac, psi_r, psi_z, B(3)
  real*8, dimension(1)                :: P, P_s, P_t, P_phi
  real*8, dimension(:), allocatable   :: P2
  real*8, dimension(:,:), allocatable :: grad_P2
  real*8  :: R_s, R_t, Z_s, Z_t, R_i, Z_i, xjac
  real*8  :: s, t, u_init_max, temp, u
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
  integer :: i_elm, i, j, k, ifail, my_id, n_cpu, ierr, n_mhd, n_geom
  real*8, dimension(fields%element_list%n_elements,2)    :: psi_minmax_list
  real*8, allocatable, dimension(:,:)             :: rans
  class(particle_base), dimension(:), allocatable :: particles_tmp
  logical, dimension(:), allocatable              :: found
  real*8  :: t0, t1
  real*8  :: Rbox(2), Zbox(2), DUMMY_R, DUMMY_Z
  integer :: blocksize, prev_blocksize, particles_to_do_local, particles_done_local
  integer :: to_find, n_tries_now, n_found
  logical :: all_done, init_uniform_space, my_include_vpar
  real*8  :: my_alpha

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

  if (present(uniform_space_rej_f)) then

    if (.not. present(uniform_space_rej_vars)) then
      write(*,*) "ERROR: if sampling function f is present variables must be given"
      call MPI_ABORT(MPI_COMM_WORLD, 10, ifail)
    end if

    ! Get the number of mhd variables to use
    allocate(P2(size(uniform_space_rej_vars,1)))

    n_mhd  = count(uniform_space_rej_vars .gt. 0)
    n_geom = size(uniform_space_rej_vars, 1) - n_mhd

    allocate(grad_P2(3,size(uniform_space_rej_vars,1)))
  
  else
    n_mhd = 0
    n_geom = 0
  end if

  if (present(alpha)) then
    write(*,*) "Not implemented yet"
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
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)

  if (.not. init_uniform_space) then
    psimin= 1d10
    psimax=-1d10
    ! Preparatory work: determine psi_min,max
    !$omp parallel do default(none) &
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
  call rng%initialize(8, random_seed(), 1, 1, ifail)

  ! Preparatory work: get R_axis, Z_axis
  call find_axis(my_id, fields%node_list, fields%element_list, psi_axis, R_axis, Z_axis, i_elm, s_axis, t_axis, ifail)

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
    ! blocks are grouped in a superblock of block*n_cpu size
    ! blocks must be equal size
    ! we communicate a blocksize here, and keep a running index of the previous starting points
    ! make this the biggest of n_tries_now
    call MPI_AllReduce(n_tries_now, blocksize, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

    call rng%jump_ahead((n_cpu-my_id-1)*prev_blocksize + my_id*blocksize)
    
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
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
#endif
    !$omp   private(i, psi, theta, phi, i_elm, s, t, R, Z, R_s, R_t, Z_s, Z_t, P2, &
    !$omp           R_i, Z_i, xjac, grad_P2, u,  &
    !$omp           P, P_s, P_t, P_phi, inv_st_jac, psi_R, psi_Z, B, H, muB, chi, ran, particle, temp, ifail, DUMMY_R, DUMMY_Z) &
    !$omp   shared(particles_tmp, psimax, psimin, found, F0, cor, mass, charge, T_Maxwell, &
    !$omp          fields, psi_minmax_list, rans, R_axis, Z_axis, blocksize, &
    !$omp          my_include_vpar, central_density, init_uniform_space, Rbox, Zbox, uniform_space_rej_vars, n_geom, n_mhd)
    do i=1,blocksize

      ran(:) = rans(:,i)

      if (init_uniform_space) then
        call transform_uniform_cylindrical([ran(3),ran(4),ran(5)], Rbox, Zbox, [0.d0,TWOPI], R, Z, phi)
        call find_RZ(fields%node_list, fields%element_list,R,Z,DUMMY_R,DUMMY_Z,i_elm,s,t,ifail)
        
        if (present(uniform_space_rej_f) .and. i_elm .ne. 0) then
    
          do k=1,n_geom
            select case (uniform_space_rej_vars(k))
            case (0);  P2(k) = 1.d0;
            case (-1); P2(k) = R;
            case (-2); P2(k) = Z;
            case (-3); P2(k) = phi;
            end select
          end do

          do k=1,n_geom
            select case (uniform_space_rej_vars(k))
            case (0);  grad_P2(:,k) = 0.d0; ! 0
            case (-1); grad_P2(:,k) = [1.d0,0.d0,0.d0]; ! R
            case (-2); grad_P2(:,k) = [0.d0,1.d0,0.d0]; ! Z
            case (-3); grad_P2(:,k) = [0.d0,0.d0,1.d0]; ! phi
            end select
          end do

          if (n_mhd .ge. 1) then
          
            call interp_PRZ(fields%node_list, fields%element_list,i_elm,                        &
                            uniform_space_rej_vars(n_geom+1:n_geom+n_mhd),n_mhd,s,t,phi,        &
                            P2(n_geom+1:n_geom+n_mhd), grad_P2(1,n_geom+1:n_geom+n_mhd),        &
                            grad_P2(2,n_geom+1:n_geom+n_mhd), grad_P2(3,n_geom+1:n_geom+n_mhd), & 
                            R_i, R_s, R_t, Z_i, Z_s, Z_t)
            
            xjac = R_s*Z_t - R_t*Z_s
            
            do k=1,n_mhd
              grad_P2(1:2,n_geom+k) = [Z_t * grad_P2(1,n_geom+k) - Z_s * grad_P2(2,n_geom+k), &
                                      -R_t * grad_P2(1,n_geom+k) + R_s * grad_P2(2,n_geom+k)]/xjac
            end do
          
          end if

          if (uniform_space_rej_f(size(uniform_space_rej_vars), P2, grad_P2) .lt. ran(7)) i_elm = 0
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

      if (i_elm .ne. 0) then
        found(i) = .true.
        ! If we are here a suitable position has been found
        particle%weight = 1. ! This is needed because the initializing values for
        ! the particles are not being used always!
        particle%i_elm = i_elm
        particle%st    = [s,t]
        particle%x     = [R,Z,phi]

        ! 1. Get B at this position
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[1],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
        inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
        psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
        psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
        ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
        B        = [+psi_Z, -psi_R, F0] / R

        ! 2. Calculate E and mu, save in particle
        call interp_PRZ(fields%node_list, fields%element_list,i_elm,[6],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)

        ! P(1)/kb/mu_zero/n_zero [K] -> multiply by kb and divide by el_chg to go to eV
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
        select type(p => particles_tmp(i))
        type is (particle_kinetic_leapfrog)

! the generic copy of particle_kinetic_leapfrog, i.e p = ..., seems broken, therefor a using yhe non-generic copy

          call copy_particle_kinetic_leapfrog( &
                 kinetic_to_kinetic_leapfrog(gc_to_kinetic(fields%node_list, fields%element_list, particle, chi, B, mass), &
                                             [0.d0, 0.d0, 0.d0], B, mass, dt=0.d0), &
                                             p )

          ! if the kinetic position is not in the grid particles(i)%i_elm the particle is lost
          if (p%i_elm .le. 0) found(i) = .false.
        type is (particle_gc)
          p = particle
        end select
      else
        found(i) = .false.
      end if
    end do
    !$omp end parallel do

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

  t_norm = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)

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
    
    if (particles(i)%i_elm .eq. 0) cycle
    
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
    particles(i)%weight = real(n * exp(-H*my_alpha/T),4) ! if zero, all particles have equal weight n
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
      particles(i)%weight = real(particles(i)%weight/P(1),4)
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
        particles(i)%weight = real(particles(i)%weight/P(1),4)
      else
        call interp_0(proj%node_list,proj%element_list,particles(i)%i_elm,[group],1, &
                      particles(i)%st(1),particles(i)%st(2),particles(i)%x(3),P)
        P = max(P, 1d-2) ! guard against divide-by-zero, maximum adjustment ratio is then 10^2
        particles(i)%weight = real(particles(i)%weight/P(1),4)
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

    if (particles(i)%i_elm .ne. 0) then
  
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


!> Dummy function to use when no transform is desired. Copies the first parameter
!> into the output or sets out to 1.
pure function no_transform(in) result(out)
  real*8, dimension(:), intent(in) :: in
  real*8 :: out
  if (size(in,1) .gt. 0) then
    out = in(1)
  else
    out = 1.d0
  end if
end function no_transform


!> Adjust weights on all particles to have the correct number of atoms in total
subroutine adjust_particle_weights(particles, num_atoms_total)
  use mpi
  class(particle_base), intent(inout), dimension(:) :: particles
  real*8, intent(in)                                :: num_atoms_total !< What the sum of the weights should be
  real*8 :: local_weights, sum_weights, local_weights_active, sum_weights_active
  integer :: ifail
 
  local_weights = sum(particles(:)%weight)

  local_weights_active = sum(particles(:)%weight, mask=particles(:)%i_elm .gt. 0)

  call MPI_AllReduce(local_weights,       sum_weights,       1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ifail)
  call MPI_AllReduce(local_weights_active,sum_weights_active,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ifail)

! Divide all weights by the sum of weights and multiply by the requested number of atoms
!  particles(:)%weight = real(real(particles(:)%weight,8) / sum_weights * num_atoms_total,4)
 
  particles(:)%weight = particles(:)%weight / sum_weights_active * num_atoms_total

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
integer :: i, ifail, seed, my_id, n_cpu
real*8, dimension(4) :: P, P_s, P_t, P_phi
real*8 :: v_out(3)
real*8 :: R, R_s, R_t, Z, Z_s, Z_t, Psi, Psi_R, Psi_Z, B(3)
real*8, parameter :: r_hat(3) = [1.d0, 0.d0, 0.d0]
real*8 :: background_kbT, background_Kelvin, background_density, V_thermal
real*8 :: DUMMY_REAL, t_norm
real*8, allocatable :: Z_coronal(:)

t_norm = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)

! Calculate a single random seed and communicate it over MPI
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ifail)
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ifail)
if (my_id .eq. 0) seed = random_seed()
call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ifail)

#ifndef WITH_Vpar
if (present(v_par) .and. v_par) then
  write(*,*) "ERROR: initialization with v// not possible with this model"
  call MPI_ABORT(-1, MPI_COMM_WORLD, ifail)
end if
#endif

do i=1,size(particles)
  if (particles(i)%i_elm .eq. 0) cycle
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
end module mod_initialise_particles
