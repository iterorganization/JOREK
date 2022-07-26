!> Tests for the CUDA kernel particle algorithms
!> These are really high level tests, there are several routines in the kernels that are not covered
!> Just the main high level routines are tested, interpolation, pushing and location
!> The tests run a threaded version of each routine and also a GPU kernel. The output of these is then
!> compared. Fields are read from an input file and particles are then initialised on this grid. The input file is
!> the same as that used for the tae_loop example.

!> The tests do not complete. This is in part due to the lack of tests of relative error in the particle_assert_equals
!> routines. Relative tests were added but have been removed as they caused the gcc compiler to crash. These tests are therefore
!> still in development. Relative checking could be done within the test routines themselves but this has not been done.
#ifdef CUDA_KERNELS
module mod_particle_kernels_test
  use fruit
  use cudafor
  use mpi
  use mod_particle_kernels
  use mod_particle_assert_equal
  use mod_particle_types,       only: particle_kinetic_leapfrog, copy_particle_kinetic_leapfrog
  use mod_particle_sim,         only: particle_group, particle_sim
  use phys_module,              only: n_particles
  implicit none

  private
  public :: run_fruit_particle_kernels

  type(particle_sim) :: sim
  integer            :: n_particles_local

contains

  subroutine run_fruit_particle_kernels()
    write(*,'(/A)') "  ... setting-up: particle kernels tests"
    call setup
    write(*,'(/A)') "  ... running: particle kernels tests"
    call test_copy_data
    call test_calc_ebpsiu_device
    call test_boris_push_cylindrical
    call test_find_rz_nearby
    call test_particle_kinetic_leapfrog_loop
    write(*,'(/A)') "  ... tearing-down: particle kernels tests"
    call teardown
  end subroutine run_fruit_particle_kernels

  !> Set-up and tear-down -------------------------
  !> set-up particle types test features
  subroutine setup()
    integer :: ierr

    !> Initialise a particle sim, we shouldn't need this but it is an easy way to set
    !> up the fields that we do need. We will ignore the particle group within the sim.
    call sim%initialize(num_groups=1)

    !> Set up the fields from a restart file. This would ideally be done in situ by calling a
    !> fields setup written for fields unti tests and not use a restart file.
    !> This approach will suffice for now but should remain under review
    call init_fields(sim)

    !> Set up particles using the H_mu_psi initialisation - this is not ideal but sufficient
    call init_particles(sim)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end subroutine setup

  subroutine teardown()

    call sim%finalize

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
    use mod_pcg32_rng

    type(particle_sim),intent(inout)  :: sim_in
    type(pcg32_rng), dimension(:), allocatable     :: rng

    real*8                                         :: rho_part, timesteps

    write(*,*) "Initialising particles using random distribution"

    ! Set up particles
    sim_in%groups(1)%Z    = 1
    sim_in%groups(1)%mass = atomic_weights(-2) !< atomic mass units

    rho_part           = 1.195d19
    timesteps          = 1e-10

    n_particles_local = int(n_particles/sim_in%n_cpu) 

    write(*,*) "Running tests with ",sim_in%n_cpu,"processes, ", n_particles ,"total particles"

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

  !> Copy data to a managed group and back and check the particles, at present the fields data is not checked
  subroutine test_copy_data

    type(particle_group_device), managed, dimension(:), allocatable   :: group_particles
    type(fields_linear_device), managed, allocatable                  :: fields
    type(particle_group), dimension(:), allocatable                   :: return_particle_groups

    write(*,*) "test copy data"

    call copy_device_data( sim , group_particles, fields )

    call copy_managed_groups_host_groups( group_particles , return_particle_groups )

    call assert_equal_particle(n_particles_local,sim%groups(1)%particles,return_particle_groups(1)%particles)

    call delete_device_data( group_particles, fields )
    deallocate( return_particle_groups(1)%particles )

    write(*,*) "test complete"

  end subroutine test_copy_data

  !> Test the interpolation with calc_EBpsiU_device
  subroutine test_calc_EBpsiU_device()
    !> variables
    real*8, dimension(8,n_particles_local) :: CPU_data, GPU_data
    real*8,parameter  :: tol_interp=7.5d-12
    integer ::  ierr

    write(*,*) sim%my_id,"test_calc_EBpsiU_device"

    !> CPU test
    call run_calc_EBpsiU_CPU(CPU_data)

    !> GPU test
    call run_calc_EBpsiU_GPU(GPU_data)

    !> Assert
    call assert_equals(CPU_data(1,:),GPU_data(1,:),n_particles_local,tol_interp,& 
        "Error calc_EBpsiU field interpolation: E direction 1 mismatch")
    call assert_equals(CPU_data(2,:),GPU_data(2,:),n_particles_local,tol_interp,& 
        "Error calc_EBpsiU field interpolation: E direction 2 mismatch")
    call assert_equals(CPU_data(3,:),GPU_data(3,:),n_particles_local,tol_interp,& 
        "Error calc_EBpsiU field interpolation: E direction 3 mismatch")
    call assert_equals(CPU_data(4,:),GPU_data(4,:),n_particles_local,tol_interp,& 
        "Error calc_EBpsiU field interpolation: B direction 1 mismatch")
    call assert_equals(CPU_data(5,:),GPU_data(5,:),n_particles_local,tol_interp,& 
        "Error calc_EBpsiU field interpolation: B direction 2 mismatch")
    call assert_equals(CPU_data(6,:),GPU_data(6,:),n_particles_local,tol_interp,& 
        "Error calc_EBpsiU field interpolation: B direction 3 mismatch")
    call assert_equals(CPU_data(7,:),GPU_data(7,:),n_particles_local,tol_interp,& 
        "Error calc_EBpsiU field interpolation: psi mismatch")
    call assert_equals(CPU_data(8,:),GPU_data(8,:),n_particles_local,tol_interp,& 
        "Error calc_EBpsiU field interpolation: U mismatch")

    write(*,*) sim%my_id,"test complete"

  end subroutine test_calc_EBpsiU_device

  subroutine run_calc_EBpsiU_CPU(data)

    !> variables
    real*8, dimension(:,:),intent(inout)     :: data

    type(particle_kinetic_leapfrog)          :: particle_tmp
    integer                                  :: i
    real*8                                   :: t, E(3), B(3), psi, U

    t = 0.0 

    select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
      do i=1,size(sim%groups(1)%particles,1)
        call copy_particle_kinetic_leapfrog(particles(i),particle_tmp)
        call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
        data(1:3,i) = E
        data(4:6,i) = B
        data(7,i)   = psi
        data(8,i)   = U
      end do
    end select

  end subroutine run_calc_EBpsiU_CPU


  subroutine run_calc_EBpsiU_GPU(data)

    real*8, dimension(:,:),intent(inout)     :: data

    type(particle_group_device), managed, dimension(:), allocatable      :: particle_groups
    type(fields_linear_device), managed, allocatable                     :: fields
    real*8, managed, dimension(:,:), allocatable                         :: data_d
    integer      :: istat, ierr_async, ierr_sync
    type(dim3)   :: grid, tBlock

    call init_gpu( sim )

    allocate(data_d(8,n_particles_local),stat=istat)
    if(istat /= 0) write(*,*) "Error allocating data"

    call copy_device_data( sim , particle_groups , fields )

    ! Launch the kernel
    tBlock = dim3(256,1,1)
    grid = dim3(ceiling(real(n_particles_local)/tBlock%x),1,1)
    call run_calc_EBpsiU_kernel<<<grid, tBlock>>>(n_particles_local, particle_groups(1), fields, data_d)
    ierr_sync = cudaGetLastError()  
    ierr_async = cudaDeviceSynchronize() 
    if (ierr_sync /= cudaSuccess) write(*,*) &
        sim%my_id, "Process. Sync kernel error:", cudaGetErrorString(ierr_sync)
    if (ierr_async /= cudaSuccess) write(*,*) &
        sim%my_id, "Process. Async kernel error:", cudaGetErrorString(ierr_async)

    !> retrieve the data
    data = data_d
    deallocate(data_d)

    call delete_device_data( particle_groups, fields )


  end subroutine run_calc_EBpsiU_GPU

  !> kernel for running calc_EBpsiU on GPUs
  attributes(global) subroutine run_calc_EBPsiU_kernel(np, group_particles, fields, data)

    integer, value, intent(in)                           :: np
    type(particle_group_device), managed, intent(inout)  :: group_particles
    type(fields_linear_device), managed , intent(inout)  :: fields
    real*8, managed, dimension(:,:),intent(inout)                 :: data

    type(particle_kinetic_leapfrog)                      :: particle_tmp
    integer                                              :: i
    real*8                                               :: t, E(3), B(3), psi, U

    i = threadIdx%x + (blockIdx%x-1) * blockDim%x
    if ( i <= np ) then
      call copy_particle_kinetic_leapfrog( group_particles%particles(i) , particle_tmp )
      if (particle_tmp%i_elm .gt. 0) then
        t = 0
        call calc_EBpsiU_device(fields, t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
        data(1,i) = E(1)
        data(2,i) = E(2)
        data(3,i) = E(3)
        data(4,i) = B(1)
        data(5,i) = B(2)
        data(6,i) = B(3)
        data(7,i) = psi
        data(8,i) = U
      endif
    endif
  end subroutine run_calc_EBPsiU_kernel

  !> Test the boris pushing algorithm
  subroutine test_boris_push_cylindrical

    type(particle_group), dimension(:), allocatable :: group_particles 
    real*8   :: start_time, particle_start_time, timestep
    integer  :: ierr
    real*8,parameter  :: tol_interp=1d-8

    write(*,*) "test_boris_push_cylindrical"

    particle_start_time = 0
    timestep = 1e-10

    call run_boris_push_cylindrical_GPU( timestep , particle_start_time , group_particles )
    call run_boris_push_cylindrical_CPU( timestep , particle_start_time )

    ! These commented lines refer to routines added to mod_particle_assert_equal that also caused gcc to crash
    ! They allow testing of relative instead of absolute values. This code was removed for this commit so as to simplify a pull request
    ! This also means that the tests as they stand fail

    !$$    call assert_equal_rel_error_particle(n_particles_local,sim%groups(1)%particles,group_particles(1)%particles)
    !$$    call set_tol(tol_interp)
    call assert_equal_particle(n_particles_local,sim%groups(1)%particles,group_particles(1)%particles)
    !$$    call reset_tol_real8()

    deallocate( group_particles(1)%particles )
    deallocate( group_particles )

    write(*,*) "test complete"

  end subroutine test_boris_push_cylindrical

  subroutine run_boris_push_cylindrical_CPU( timestep, particle_start_time )

    use mpi
    use omp_lib
    use mod_boris, only: boris_push_cylindrical

    real*8, intent(in)                       :: timestep, particle_start_time

    !> variables
    type(particle_kinetic_leapfrog)          :: particle_tmp
    integer                                  :: nthreads, i
    real*8                                   :: t, E(3), B(3), psi, U

    t = particle_start_time

    if( nthreads <= 0 .or. nthreads > omp_get_max_threads()) nthreads = omp_get_max_threads()

    select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
      !$omp parallel do default(shared) &
      !$omp private(particle_tmp,i,E,B,psi,U) &
      !$omp num_threads(nthreads) &
      !$omp schedule(dynamic,10)
      do i=1,size(sim%groups(1)%particles,1)
        call copy_particle_kinetic_leapfrog(particles(i),particle_tmp)
        if (particle_tmp%i_elm .gt. 0) then
          call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
          call boris_push_cylindrical(particle_tmp, sim%groups(1)%mass, E, B, timestep)
        endif
        call copy_particle_kinetic_leapfrog(particle_tmp,particles(i))
      enddo
      !$omp end parallel do
    end select

  end subroutine run_boris_push_cylindrical_CPU

  subroutine run_boris_push_cylindrical_GPU( timestep, particle_start_time , return_particle_groups )

    real*8, intent(inout)                                                :: particle_start_time, timestep
    type(particle_group), dimension(:), allocatable, intent(inout)       :: return_particle_groups

    type(particle_group_device), managed, dimension(:), allocatable      :: particle_groups
    type(fields_linear_device), managed, allocatable                     :: fields
    integer      :: istat
    type(dim3)   :: grid, tBlock

    call init_gpu( sim )

    ! Copy the data
    call copy_device_data( sim , particle_groups , fields )

    ! Launch the kernel
    tBlock = dim3(256,1,1)
    grid = dim3(ceiling(real(n_particles_local)/tBlock%x),1,1)
    call run_boris_push_cylindrical_kernel<<<grid, tBlock>>>(n_particles_local, particle_groups(1), fields, timestep , particle_start_time)
    istat = cudaDeviceSynchronize()
    if (istat /= cudaSuccess) write(*,*) cudaGetErrorString(istat)
    call copy_managed_groups_host_groups( particle_groups, return_particle_groups )

    call delete_device_data( particle_groups, fields )

  end subroutine run_boris_push_cylindrical_GPU

  !> kernel for running calc_EBpsiU on GPUs
  attributes(global) subroutine run_boris_push_cylindrical_kernel(np, group_particles, fields, timestep, particle_start_time)

    use mod_boris, only: boris_push_cylindrical

    integer, value, intent(in)                           :: np
    type(particle_group_device), managed, intent(inout)  :: group_particles
    type(fields_linear_device), managed , intent(inout)  :: fields
    real*8, value, intent(in)                            :: timestep, particle_start_time

    type(particle_kinetic_leapfrog)                      :: particle_tmp      
    integer                                              :: i
    real*8                                               :: t, E(3), B(3), psi, U

    t = particle_start_time
    i = threadIdx%x + (blockIdx%x-1) * blockDim%x 
    if ( i <= np ) then
      call copy_particle_kinetic_leapfrog( group_particles%particles(i) , particle_tmp )
      if (particle_tmp%i_elm .gt. 0) then
        call calc_EBpsiU_device(fields, t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
        call boris_push_cylindrical(particle_tmp, group_particles%mass, E, B, timestep)
      endif
      call copy_particle_kinetic_leapfrog( particle_tmp, group_particles%particles(i) )
    endif

  end subroutine run_boris_push_cylindrical_kernel

  !> Test find of element index
  subroutine test_find_rz_nearby()
    !> variables
    integer, dimension(n_particles_local) :: CPU_data, GPU_data
    integer :: i

    write(*,*) "test_find_rz_nearby"

    ! Add 1 to i_elm so that find_rz has to find the particle again
    do i = 1, n_particles_local
      sim%groups(1)%particles(i)%i_elm = sim%groups(1)%particles(i)%i_elm + 1
    end do

    !> GPU test - run this first as the original sim type is unaffected
    call run_find_rz_nearby_GPU(GPU_data)

    !> CPU test
    call run_find_rz_nearby_CPU(CPU_data)

    !> Assert
    call assert_equals(CPU_data(:),GPU_data(:),n_particles_local,& 
        "Error find_rz_nearby: element index mismatch")

    do i = 1, n_particles_local
      sim%groups(1)%particles(i)%i_elm = sim%groups(1)%particles(i)%i_elm - 1
    end do

    write(*,*) "test complete"

  end subroutine test_find_rz_nearby

  subroutine run_find_rz_nearby_CPU(data)
    use mod_find_rz_nearby, only: find_rz_nearby
    !> variables
    integer, dimension(:),intent(inout)     :: data

    type(particle_kinetic_leapfrog)          :: particle_tmp
    integer                                  :: i, i_elm_old, ifail, i_elm_new
    real*8                                   :: rz_old(2), st_old(2)

    select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
      do i=1,size(sim%groups(1)%particles,1)
        call copy_particle_kinetic_leapfrog(particles(i),particle_tmp)
        rz_old    = particle_tmp%x(1:2)
        st_old    = particle_tmp%st
        i_elm_old = particle_tmp%i_elm
        call find_rz_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), &
            st_old(1), st_old(2), i_elm_old, particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), &
            particle_tmp%st(2), i_elm_new, ifail)
        data(i) = i_elm_new
      end do
    end select

  end subroutine run_find_rz_nearby_CPU


  subroutine run_find_rz_nearby_GPU(data)

    integer, dimension(:),intent(inout)     :: data

    type(particle_group_device), managed, dimension(:), allocatable   :: particle_groups
    type(fields_linear_device), managed, allocatable                  :: fields
    integer, managed, dimension(:), allocatable                       :: data_d      
    integer      :: istat, ierrSync, ierrAsync
    type(dim3)   :: grid, tBlock

    ! Set the device to use for this process 
    call init_gpu( sim )

    allocate(data_d(n_particles_local))

    ! Copy the data
    call copy_device_data( sim , particle_groups , fields )

    ! Launch the kernel
    tBlock = dim3(256,1,1)
    grid = dim3(ceiling(real(n_particles_local)/tBlock%x),1,1)
    call run_find_rz_nearby_kernel<<<grid, tBlock>>>(n_particles_local, particle_groups(1), fields, data_d)
    ierrSync = cudaGetLastError()
    ierrAsync = cudaDeviceSynchronize()
    if (ierrSync /= cudaSuccess) write(*,*) &
        "Sync kernel error:", cudaGetErrorString(ierrSync)
    if (ierrAsync /= cudaSuccess) write(*,*) &
        "Async kernel error:", cudaGetErrorString(ierrAsync)

    !> retrieve the data
    data = data_d
    deallocate(data_d)

    call delete_device_data( particle_groups, fields )

  end subroutine run_find_rz_nearby_GPU

  !> kernel for running calc_EBpsiU on GPUs
  attributes(global) subroutine run_find_rz_nearby_kernel(np, group_particles, fields, data)
    use mod_find_rz, only: find_rz_nearby_device
    integer, value, intent(in)                           :: np
    type(particle_group_device), managed, intent(inout)  :: group_particles
    type(fields_linear_device), managed , intent(inout)  :: fields
    integer, managed, dimension(:),intent(inout)         :: data

    type(particle_kinetic_leapfrog)                      :: particle_tmp      
    integer                                              :: i, i_elm_old, i_elm_new, ifail
    real*8                                               :: rz_old(2), st_old(2)

    i = threadIdx%x + (blockIdx%x-1) * blockDim%x 
    if ( i <= np ) then
      call copy_particle_kinetic_leapfrog( group_particles%particles(i) , particle_tmp )
      if(particle_tmp%i_elm .gt. 0) then
        rz_old    = particle_tmp%x(1:2)
        st_old    = particle_tmp%st
        i_elm_old = particle_tmp%i_elm
        call find_rz_nearby_device(fields%node_list,fields%element_list,rz_old(1),rz_old(2), &
            st_old(1),st_old(2),i_elm_old,particle_tmp%x(1), particle_tmp%x(2), &
            particle_tmp%st(1),  particle_tmp%st(2), i_elm_new, ifail)
        data(i) = i_elm_new       
      end if
    end if

  end subroutine run_find_rz_nearby_kernel

  !> Test the full particle loop
  subroutine test_particle_kinetic_leapfrog_loop

    type(particle_group), dimension(:), allocatable :: group_particles 
    real*8   :: start_time, particle_start_time, timestep, tol
    integer  :: n_steps, ierr

    tol = 5e-10

    write(*,*) "test_particle_kinetic_leapfrog_loop"

    particle_start_time = 0
    n_steps = 1000
    timestep = 1e-10

    start_time = MPI_Wtime()
    call particle_kinetic_leapfrog_full_loop( sim , n_steps , timestep , particle_start_time , group_particles)
    write(*,*) "Proc",sim%my_id, "GPU full loop completed in ",MPI_Wtime()-start_time," s"

    start_time = MPI_Wtime()
    call run_particle_kinetic_leapfrog_loop_CPU( n_steps, timestep, particle_start_time )
    write(*,*) "Proc",sim%my_id, "Threaded full loop completed in ",MPI_Wtime()-start_time," s"

    ! These commented lines refer to routines added to mod_particle_assert_equal that also caused gcc to crash
    ! They allow testing of relative instead of absolute values. This code was removed for this commit so as to simplify a pull request
    ! This also means that the tests as they stand fail
!!$       call set_tol(tol)
!!$       call assert_equal_rel_error_particle(n_particles_local,sim%groups(1)%particles,group_particles(1)%particles)
!!$       call reset_tol_real8()

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    deallocate( group_particles(1)%particles )
    deallocate( group_particles )

    write(*,*) "test complete"

  end subroutine test_particle_kinetic_leapfrog_loop

  subroutine run_particle_kinetic_leapfrog_loop_CPU(n_steps, timestep, particle_start_time )
    use phys_module, only: F0, tstep
    use mpi
    use omp_lib
    use mod_find_rz_nearby, only: find_rz_nearby
    use mod_boris, only: boris_push_cylindrical

    integer, intent(in)                      :: n_steps
    real*8, intent(in)                       :: timestep, particle_start_time

    !> variables
    type(particle_kinetic_leapfrog)          :: particle_tmp
    integer                                  :: i, j, ifail, i_elm_old
    real*8                                   :: t, rz_old(2), st_old(2), E(3), B(3), psi, U, start_time

    select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
      !$omp parallel default(shared) &
      !$omp private(particle_tmp,i,j,E,B,psi,U,rz_old,st_old,i_elm_old,ifail) 
      !$omp do ! schedule(dynamic,10)
      do i=1,n_particles_local
        call copy_particle_kinetic_leapfrog(particles(i),particle_tmp)            
        do j=1,n_steps
          if (particle_tmp%i_elm .le. 0) exit             
          t = particle_start_time + (j-1)*timestep
          call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
          call boris_push_cylindrical(particle_tmp, sim%groups(1)%mass, E, B, timestep)               
          if (particle_tmp%i_elm .gt. 0) then
            rz_old    = particle_tmp%x(1:2)
            st_old    = particle_tmp%st
            i_elm_old = particle_tmp%i_elm                  
            call find_rz_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), &
                st_old(1), st_old(2), i_elm_old, particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), &
                particle_tmp%st(2), particle_tmp%i_elm, ifail)                  
          endif
        end do
        call copy_particle_kinetic_leapfrog(particle_tmp,particles(i))
      enddo
      !$omp end do
      !$omp end parallel
    end select

  end subroutine run_particle_kinetic_leapfrog_loop_CPU


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
