#ifdef CUDA_KERNELS
!> Contains routines for offloading particle computation to GPU devices using CUDA kernels
!> At present only the particle_kinetic_leapfrog type is supported

module mod_particle_kernels

  use cudafor
  use mpi
  use mod_particle_sim, only: particle_sim
  use mod_particle_types, only: particle_kinetic_leapfrog, copy_particle_kinetic_leapfrog
  use data_structure
  
  implicit none

  type fields_linear_device
     type(type_node_list)       :: node_list        !< Current node list
     type(type_element_list)    :: element_list     !< Current element list
     logical                    :: static           !< if true do not time interpolate
     logical                    :: flag_zero_dpsidt !< if true, P_time(1) = dpsi/dt = 0
     real*8                     :: time_now         !< Time of current restart file (SI units)
     real*8                     :: time_prev        !< Time of previous restart file (SI units)
  end type fields_linear_device

  private
  public fields_linear_device, particle_kinetic_leapfrog_loop, copy_device_data
  
contains


  subroutine particle_kinetic_leapfrog_loop( sim )! , n_steps , timestep , particle_start_time )
    type(particle_sim), intent(inout)                                    :: sim
    
    type(particle_kinetic_leapfrog), managed, dimension(:), allocatable  :: particles
    type(fields_linear_device), managed                                         :: fields
    integer      :: n_particles, tBlock_size, istat, n_steps
    type(dim3)   :: grid, tBlock

    ! Print out some info on the GPU
    if(sim%my_id==0) call device_query()

    n_steps = 1000
    n_particles = size( sim%groups(1)%particles,1)
    tBlock_size = 256
    tBlock = dim3(tBlock_size,1,1)
    grid = dim3(n_particles/tBlock_size,1,1)

    istat = cudaSetDevice(sim%my_id)
    if (istat /= cudaSuccess) write(*,*) cudaGetErrorString(istat)

    call copy_device_data( sim , particles , fields )

    call particle_kinetic_leapfrog_loop_kernel<<<grid, tBlock>>>(n_particles, particles, fields, n_steps) ! timesteps, particle_start_time )
    
  end subroutine particle_kinetic_leapfrog_loop

  attributes(global) subroutine particle_kinetic_leapfrog_loop_kernel(n_particles, particles, fields, n_steps) !, nodes, elements, fields_meta, n_steps, timesteps, particle_start_time )
    type(particle_kinetic_leapfrog), managed, dimension(:), intent(inout)  :: particles
    type(fields_linear_device), managed , intent(inout)                           :: fields
    integer, value, intent(in)                                             :: n_particles, n_steps

    type(particle_kinetic_leapfrog)         :: particle_tmp
    integer                                 :: i,j
   
    i = threadIdx%x + (blockIdx%x-1) * blockDim%x 
    if(i == 10) write(*,*) "Seems to be working"
    if ( i <= n_particles ) then
       call copy_particle_kinetic_leapfrog( particles(i) , particle_tmp )
       do j=1,n_steps
          if (particle_tmp%i_elm .le. 0) then
             write(*,*) "Losing particle",j
             exit
          endif
       end do
    end if
    
  end subroutine particle_kinetic_leapfrog_loop_kernel
  
  subroutine device_query()

    type(cudaDeviceProp) :: prop    
    integer              :: istat, i, n_devices
    
    istat = cudaGetDeviceCount(n_devices)
    do i = 0, n_devices-1
       istat = cudaGetDeviceProperties(prop, i)
       write(*,"(' Device Number: ',i0)") i
       write(*,"('   Device name: ',a)") trim(prop%name)
       write(*,"('   Compute capability: ', i0,'.', i0)") prop%major,prop%minor
       write(*,"('   Max threads per block: ', i0)") prop%maxThreadsPerBlock
       write(*,"('   Number of multiprocessors: ', i0)") prop%multiProcessorCount
       write(*,"('   Number threads per multiprocessor: ', i0)") prop%maxThreadsPerMultiProcessor
       write(*,"('   Memory Clock Rate (KHz): ', i0)") &
            prop%memoryClockRate
       write(*,"('   Memory Bus Width (bits): ', i0)") &
            prop%memoryBusWidth
       write(*,"('   Peak Memory Bandwidth (GB/s): ', f6.2)") &
            2.0*prop%memoryClockRate*(prop%memoryBusWidth/8)/10.0**6
       write(*,*)        
    enddo

  end subroutine device_query

  subroutine copy_device_data( sim , particles , fields )

    type(particle_sim), intent(inout)                                                  :: sim
    type(particle_kinetic_leapfrog), managed, dimension(:), allocatable, intent(inout) :: particles
    type(fields_linear_device), managed , intent(inout)                                       :: fields
    
    write(*,*) "Copying device data"

    call copy_particles( sim , particles )
    
    call copy_fields_device( sim , fields )
    
  end subroutine copy_device_data


  subroutine copy_particles(sim , particles)    
    type(particle_sim), intent(inout)                                                  :: sim
    type(particle_kinetic_leapfrog), managed, dimension(:), allocatable, intent(inout) :: particles

    integer      :: i, n_particles

    n_particles = size( sim%groups(1)%particles,1)

    allocate(particles(n_particles))
    
    select type(p => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
       do i=1,n_particles
          call copy_particle_kinetic_leapfrog( p(i) , particles(i) )
       end do
    end select

  end subroutine copy_particles
  
  subroutine copy_fields_device( sim , fields )
    type(particle_sim), intent(inout)    :: sim
    type(fields_linear_device), intent(inout)   :: fields
    
    call copy_element_list( sim%fields%element_list, fields%element_list )

  end subroutine copy_fields_device

  subroutine copy_element_list( in , out )
    type(type_element_list), intent(in)        :: in
    type(type_element_list), intent(inout)     :: out

    integer  :: i
    
    do i=1,n_elements_max
       call copy_element( in%element(i) , out%element(i) )
    end do
   
  end subroutine copy_element_list

  !> Copy and element type from the host to a GPU device
  subroutine copy_element( in , out ) 

    type(type_element), intent(in)  :: in 
    type(type_element), intent(out) :: out 

    out%vertex         = in%vertex
    out%neighbours     = in%neighbours
    out%size           = in%size  
    out%father         = in%father
    out%n_sons         = in%n_sons  
    out%n_gen          = in%n_gen 
    out%sons           = in%sons
    out%contain_node   = in%contain_node
    out%nref           = in%nref
    
  end subroutine copy_element

end module mod_particle_kernels
#endif
