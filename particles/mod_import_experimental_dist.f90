module mod_import_experimental_dist
  use mod_rng
  use particle_tracer
  use data_structure
  use mod_particle_types
  use constants
  use mod_interp
  use iso_c_binding
  use hdf5_io_module
  use hdf5
  implicit none
  private
  public import_particles, calculate_B

contains

! > Subroutine to import particles from a 4D distribution function:
! >> particles:             list of particles from simulation
! >> fields:                JOREK fields on which to initialise particles
! >> filename:              hdf5 filename containing 4D distribution function
! >> mass:                  mass in amu of particles
! >> n_phi_planes:          number of times to copy particle with increasing phi (2 pi / n_phi_planes) increment,
!                            should be chosen 2x highest toroidal harmonic if used to prevent aliasing.
! >> fraction_phi_planes:   fraction of particles that is copied with increasing phi, e.g. 0.99d0 will ensure 99% of 
!                            the particles are copied with increasing phi (and thus provide no toroidal harmonics).
!                            Noise in harmonics follows n0*sqrt(1-fraction_phi_planes) dependence (n0 noise with this 0.d0)
! For documentation of the HDF5 file: https://www.jorek.eu/wiki/doku.php?id=experimental_fast_particle
subroutine import_particles(particles,fields,filename, rng_base,mass, n_phi_planes_in,fraction_phi_planes)
  use mpi
  use mod_sampling
  use mod_random_seed
  use mod_interp
!$ use omp_lib

  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  class(fields_base),    intent(in)                 :: fields
  class(type_rng),       intent(in)                 :: rng_base !< What type of random number generator to use (will be reseeded here)
  character*(*),         intent(in)                 :: filename
  real*8,                intent(in)                 :: mass, fraction_phi_planes
  integer,               intent(in)                 :: n_phi_planes_in

  integer                                           :: my_id, n_cpu, ierr, rank, n_particles,ifail
  integer                                           :: hdferr
  integer, allocatable, dimension(:)                :: particles_per_proc
  integer*8                                         :: arraysize(4),maxdims(4)
  integer*8                                         :: E1Dsize(1),R1Dsize(1),Z1Dsize(1),Pitch1Dsize(1),tmp(1)
  integer*8                                         :: OneD_arraysize
  real*8,  allocatable, dimension(:,:,:,:)          :: F0_norm
  real*8, allocatable, dimension(:)                 :: E1D,Pitch1D,R1D,Z1D

  character(len=1024) :: myidfile

  integer(HID_T)    :: file,dset, file_space

  ! Particles
  class(type_rng), allocatable                      :: rng
  type(particle_gc)                                 :: particle
  integer                                           :: particles_to_do_local, particles_done_local, prev_blocksize,to_find, n_phi_planes_dummy
  integer                                           :: n_tries_now, blocksize,i, i_elm
  integer                                           :: index_R,index_Z,index_E,index_Pitch, n_found, i_found, i_particle,i_plane
  logical                                           :: all_done
  logical, dimension(:), allocatable                :: found
  real*8, allocatable, dimension(:,:)               :: rans
  real*8                                            :: Rbox(2), Zbox(2),R,Z,phi, ran(8), dummy_R, dummy_Z, PitchBox(2), EBox(2)
  real*8                                            :: s,t,Pitch_tmp, E_tmp, grids_start(4), grids_spacing(4), F0_value, mu_tmp
  real*8                                            :: B_tmp(3)
  class(particle_base), dimension(:), allocatable   :: particles_tmp

  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)
  if(my_id .eq. 0 )then
    write(*,*) "*************************************"
    write(*,*) "* Initialising from experimental    *"
    write(*,*) "* fast ion distribution function    *"
    write(*,*)"*************************************"
    write(*,*)" "
  endif
  n_phi_planes_dummy = n_phi_planes_in
  ! First read the HDF5 with one MPI task, broadcast to all other

  if(my_id .eq. 0 ) then
    call h5open_f(hdferr)
    write(*,*) "PARTICLES: Opening '", filename, "' for initialising particles"
    write(*,*) "PARTICLES: from realistic distribution function"
    call h5fopen_f(filename,h5f_acc_rdonly_f,file,hdferr)
    call h5dopen_f(file, "F0_norm", dset, hdferr)
    call h5dget_space_f(dset, file_space, hdferr)

    call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
    ! 2 spatial (R,Z), 2 velocity space (Pitch angle & E). Gyrophase and phi are taken to be uniformly random.
    if(rank .ne. 4 ) then
      write(*,*) "PARTICLES: rank of realistic distribution function should be 4, while it is ", rank,"."
      write(*,*) "PARTICLES: Aborting..."
      call MPI_ABORT(MPI_COMM_WORLD, -1,ierr)
    endif

    call h5sget_simple_extent_dims_f(file_space, arraysize, maxdims, hdferr)
    write(*,"(A,4I4)") " PARTICLES: Dimensionality of F0 is ", arraysize

    call h5sclose_f(file_space, hdferr)
    call h5dclose_f(dset, hdferr)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_BCAST(arraysize,4,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
  allocate(F0_norm(arraysize(1),arraysize(2),arraysize(3),arraysize(4)))

  if (my_id .eq. 0) then
    ! Read F0_norm
    call HDF5_array4D_reading(file,F0_norm,"F0_norm")

    call h5dopen_f(file, "R_1D", dset, hdferr)
    call h5dget_space_f(dset, file_space, hdferr)
    call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
    call h5sget_simple_extent_dims_f(file_space, R1Dsize, tmp, hdferr)
    call h5sclose_f(file_space, hdferr)
    call h5dclose_f(dset, hdferr)
    write(*,"(A,1I4)") " PARTICLES: R1D size     = ", R1Dsize

    call h5dopen_f(file, "Z_1D", dset, hdferr)
    call h5dget_space_f(dset, file_space, hdferr)
    call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
    call h5sget_simple_extent_dims_f(file_space, Z1Dsize, tmp, hdferr)
    call h5sclose_f(file_space, hdferr)
    call h5dclose_f(dset, hdferr)
    write(*,"(A,1I4)") " PARTICLES: Z1D size     = ", Z1Dsize

    call h5dopen_f(file, "Pitch_1D", dset, hdferr)
    call h5dget_space_f(dset, file_space, hdferr)
    call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
    call h5sget_simple_extent_dims_f(file_space, Pitch1Dsize, tmp, hdferr)
    call h5sclose_f(file_space, hdferr)
    call h5dclose_f(dset, hdferr)
    write(*,"(A,1I4)") " PARTICLES: Pitch1D size = ", Pitch1Dsize

    call h5dopen_f(file, "E_1D", dset, hdferr)
    call h5dget_space_f(dset, file_space, hdferr)
    call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
    call h5sget_simple_extent_dims_f(file_space, E1Dsize, tmp, hdferr)
    call h5sclose_f(file_space, hdferr)
    call h5dclose_f(dset, hdferr)
    write(*,"(A,1I4)") " PARTICLES: E1D size     = ", E1Dsize
  endif
  call MPI_BCAST(R1Dsize,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(E1Dsize,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Z1Dsize,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Pitch1Dsize,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
  allocate(R1D(R1Dsize(1)))
  allocate(Z1D(Z1Dsize(1)))
  allocate(Pitch1D(Pitch1Dsize(1)))
  allocate(E1D(E1Dsize(1)))
  if(my_id .eq. 0) then

    call HDF5_array1D_reading(file,R1D,"R_1D")
    call HDF5_array1D_reading(file,Z1D,"Z_1D")
    call HDF5_array1D_reading(file,Pitch1D,"Pitch_1D")
    call HDF5_array1D_reading(file,E1D,"E_1D")

    ! Assuming the array in Python is done by R,Z,Pitch,Energy, this will be reverse here. Sanity check:

    if ( R1Dsize(1) .ne. arraysize(4) .or. Z1Dsize(1) .ne. arraysize(3) .or. Pitch1Dsize(1) .ne. arraysize(2) .or. E1Dsize(1) .ne. arraysize(1)) then
      write(*,*) "PARTICLES: Dimensions of 1D arrays and F0 do not conform."
      write(*,*) "PARTICLES: Please check correctness of H5 file. Aborting..."
      write(*,*) "PARTICLES: R    ", R1Dsize, arraysize(4)
      write(*,*) "PARTICLES: Z    ", Z1Dsize, arraysize(3)
      write(*,*) "PARTICLES: Pitch", Pitch1Dsize, arraysize(2)
      write(*,*) "PARTICLES: E    ", E1Dsize, arraysize(1)
      call MPI_ABORT(MPI_COMM_WORLD,10,ierr)
    endif
    call h5fclose_f(file,hdferr)
    call h5close_f(hdferr)
    write(*,*) "PARTICLES: Succesfully read in the realistic distribution function"

  endif

  write(*,"(4I4)"), R1Dsize,Z1Dsize,E1Dsize,Pitch1Dsize
  OneD_arraysize = Z1Dsize(1)*Pitch1Dsize(1)*E1Dsize(1)*R1Dsize(1)
  write(*,"(A,I9,A,I3)") " PARTICLES: total F0 elements = ", OneD_arraysize, " on process ",my_id
  ! Starting broadcasting....
  call MPI_BCAST(F0_norm,int(OneD_arraysize,4), MPI_REAL8, 0 , MPI_COMM_WORLD, ierr)
  call MPI_BCAST(R1D,int(R1Dsize(1),4),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Z1D,int(Z1Dsize(1),4),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(E1D,int(E1Dsize(1),4),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Pitch1D,int(Pitch1Dsize(1),4),MPI_REAL8,0,MPI_COMM_WORLD,ierr)



  ! Particles.
  write(*,"(A,I9,A,I3)") " PARTICLES: initialising ", size(particles)," on process ", my_id
  allocate(rng,source=rng_base)
  call rng%initialize(8, random_seed(), 1, 1, ifail)
  particles(:)%i_elm = 0

  particles_to_do_local  = size(particles)
  particles_done_local   = 0
  prev_blocksize = 0 ! initial block size
  all_done = .false.

  ! Boxes for rejection sampling
  call domain_bounding_box(fields%node_list, fields%element_list, Rbox(1), Rbox(2), Zbox(1), Zbox(2))
  PitchBox(1)=minval(Pitch1D)
  PitchBox(2)=maxval(Pitch1D)
  EBox(1)=minval(E1D)
  EBox(2)=maxval(E1D)

  grids_start = [R1D(1),Z1D(1),Pitch1D(1),E1D(1)]
  grids_start = [E1D(1),Pitch1D(1),Z1D(1),R1D(1)]
  grids_spacing = [E1D(2)-E1D(1),Pitch1D(2)-Pitch1D(1),Z1D(2)-Z1D(1),R1D(2)-R1D(1)]


  do while (.not. all_done)
    to_find = (particles_to_do_local-particles_done_local)
    n_tries_now = to_find
    if (n_tries_now .lt. 64 .and. n_tries_now .gt. 0) n_tries_now = 64

    call MPI_AllReduce(n_tries_now, blocksize, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call rng%jump_ahead((n_cpu-my_id-1)*prev_blocksize + my_id*blocksize)
    prev_blocksize = blocksize
    allocate(rans(8,blocksize))
    do i=1,blocksize
      call rng%next(rans(:,i))
    end do

    ! Allocate temporary particle storage.
    select type (particles)
      type is (particle_kinetic)
        write(*,*) "ERROR: particle_kinetic not supported yet for import_particles"
      type is (particle_kinetic_leapfrog)
        allocate(particle_kinetic_leapfrog::particles_tmp(blocksize))
      type is (particle_gc)
        allocate(particle_gc::particles_tmp(blocksize))
      class default
        write(*,*) "ERROR: particle type not supported yet for import_particles"
        call exit(1)
    end select
    allocate(found(blocksize))

    ! OMP loop for particles
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
    !$omp shared(blocksize,rans,Rbox,Zbox,fields,F0_norm,PitchBox, EBox,my_id,grids_start,grids_spacing,found,R1D,Z1D,&
    !$omp        Z1Dsize,R1Dsize,particles_tmp,mass)&
#endif
    !$omp private(i,R,Z,phi,ran,dummy_R,dummy_Z,i_elm,s,t,ifail,Pitch_tmp, E_tmp,ierr,F0_value,particle,B_tmp)
    do i=1,blocksize
      ran(:) = rans(:,i)
      found(i) = .false.
      R = (Rbox(2)-Rbox(1))*ran(1)+Rbox(1)
      Z = (Zbox(2)-Zbox(1))*ran(2)+Zbox(1)
      !call find_RZ(fields%node_list, fields%element_list,R,Z,dummy_R,dummy_Z,i_elm,s,t,ifail)
      ! At  this point, we have R,Z and a check if the particle is in the domain.
      ! Now, we need to generate numbers corresponding to pitch angle and energy.
      Pitch_tmp = (PitchBox(2)-PitchBox(1))*ran(4)+PitchBox(1)
      E_tmp = (EBox(2)-EBox(1))*ran(5)+EBox(1)

      ! We have now 4 numbers corresponding to R,Z, Pitch Angle. Do a 4D linear interpolation.
      if(R > R1D(1).and. Z> Z1D(1) .and. R < R1D(R1Dsize(1)) .and. Z < Z1D(Z1Dsize(1))) then
        F0_value = FourD_linear_interp(F0_norm,grids_start,grids_spacing,[E_tmp,Pitch_tmp,Z,R])
        ! Accept particle if ran(5) < F0_value (4D rejection sampled)
        if(ran(6)< F0_value) then


          call find_RZ(fields%node_list, fields%element_list,R,Z,DUMMY_R,DUMMY_Z,i_elm,s,t,ifail)
          if(i_elm > 0) then
            found(i) = .true.
          else
            cycle
          endif
          particle%x = [R,Z, ran(3)*TWOPI]
          particle%weight = 1
          particle%i_elm  = i_elm
          particle%st     = [ s, t ]
          particle%E= E_tmp ! in eV

          B_tmp = calculate_B(fields,i_elm,s,t,particle%x(3))
          particle%mu = sign(E_tmp*(1- Pitch_tmp**2)/norm2(B_tmp),-1.d0*Pitch_tmp) ! Pitch is wrt current here, current opposite to B.

          particle%q     = int(1,1)

          select type(p => particles_tmp(i))
            type is (particle_kinetic_leapfrog)
              call copy_particle_kinetic_leapfrog( &
              kinetic_to_kinetic_leapfrog(gc_to_kinetic(fields%node_list, fields%element_list, particle, TWOPI*ran(7), B_tmp, mass), &
              [0.d0, 0.d0, 0.d0], B_tmp, mass, dt=0.d0), p )
              if(p%i_elm < 1) found(i)=.false.
            type is(particle_gc)
              p = particle
          end select
        endif

      endif ! i_elm > 0
    enddo
    !$omp end parallel do
    n_found = count(found)



    i_found=1
    do i_particle=1,size(particles_tmp)
      if(found(i_particle)) then
        if(real(particles_done_local+i_found)/real(particles_to_do_local)>fraction_phi_planes .and. n_phi_planes_dummy .ne. 1) then
          n_phi_planes_dummy = 1
          write(*,"(A,F8.1,A,I4)") " PARTICLES: Set n_phi_planes = 1 at", real(particles_done_local+i_found)/real(particles_to_do_local)*100.d0,"% at process", my_id
        endif
        do i_plane = 1,n_phi_planes_dummy
          if(i_found .gt. to_find) exit
          select type(p1 => particles(particles_done_local+i_found))
            type is(particle_kinetic_leapfrog)
              select type(p2 => particles_tmp(i_particle))
                type is (particle_kinetic_leapfrog)
                  p2%x(3) = p2%x(3)+TWOPI/real(n_phi_planes_in,8)
                  p1 = p2
              end select
            type is(particle_gc)
              select type(p2 => particles_tmp(i_particle))
                type is (particle_gc)
                  p2%x(3) = p2%x(3)+TWOPI/real(n_phi_planes_in,8)
                  p1 = p2
              end select
          end select
          i_found=i_found+1
        enddo
      endif

    enddo

    particles_done_local = particles_done_local + i_found-1

    call MPI_AllReduce(particles_to_do_local .eq. particles_done_local, all_done, &
      1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    deallocate(rans,found,particles_tmp)
    write(*,"(A,I3,A,I10,A,I10,A,F5.1,A)") " PARTICLES: process", my_id, " found ", n_found, " out of", to_find,". At ", real(particles_done_local)/real(particles_to_do_local)*100.d0,"%."
  enddo
  ! De-allocating F0 and grids on each MPI process.

  deallocate(F0_norm)
  deallocate(R1D)
  deallocate(Z1D)
  deallocate(E1D)
  deallocate(Pitch1D)

end subroutine import_particles

pure function FourD_linear_interp(F0_norm,starts,spacings,values) result(value_out)
  real*8, dimension(:,:,:,:), intent(in) :: F0_norm
  real*8,                     intent(in) :: values(4),starts(4),spacings(4)
  real*8                                 :: up_fractions(4), down_fractions(4)
  integer                                :: minimum_indices(4)
  real*8                                 :: value_out
  integer                                :: i,j,k,l
  integer                                :: index_tmp(4)
  minimum_indices = floor((values - starts)/spacings)+1
  ! Weight fraction to lower gridpoint (i.e. [0,1] with value 0.3 will yield 0.7 weight)
  up_fractions = mod(values-starts,spacings)/spacings
  down_fractions = 1.d0 - up_fractions
  value_out = 0.d0
  ! Visualize as the corners of a 4D hypercube. Then just go along each of the 16 points and add up each
  ! contribution.
  do i=0,1
    do j=0,1
      do k=0,1
        do l = 0,1
          index_tmp = minimum_indices+[i,j,k,l]
          value_out = value_out +F0_norm(index_tmp(1),index_tmp(2),index_tmp(3),index_tmp(4))*product(real([i,j,k,l],8)*up_fractions+real(1-[i,j,k,l],8)*down_fractions)
        enddo
      enddo
    enddo
  enddo

end function

pure function calculate_B(fields, i_elm,s,t,phi) result(B)
  use phys_module, only: F0
  class(fields_base),                  intent(in)    :: fields
  integer,                             intent(in)    :: i_elm
  real*8,                              intent(in)    :: s,t,phi

#ifdef fullmhd
  real*8                    :: A3, AR, AZ, A3_R, A3_Z, AR_Z, AR_p, AZ_R, AZ_P, Fprof
  real*8, dimension(3)      :: P, P_s, P_t, P_phi
#else
  real*8, dimension(1)      :: P, P_s, P_t, P_phi
#endif
  real*8                                             :: B(3)
  real*8                                             :: R, Z, inv_st_jac, psi_r, psi_z
  real*8                                             :: R_s, R_t, Z_s, Z_t, R_i, Z_i, xjac

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

#else
  call interp_PRZ(fields%node_list, fields%element_list,i_elm,[1],1,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)
  inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
  psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
  psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
  ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
  B        = [+psi_Z, -psi_R, F0] / R
#endif

end function

end module mod_import_experimental_dist


