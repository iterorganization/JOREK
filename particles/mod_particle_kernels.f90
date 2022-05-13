#ifdef CUDA_KERNELS
!> Contains routines for offloading particle computation to GPU devices using CUDA kernels
!> At present only the particle_kinetic_leapfrog type is supported

module mod_particle_kernels
  use cudafor
  use mpi
  use mod_particle_sim,   only: particle_sim, particle_group
  use mod_particle_types, only: particle_kinetic_leapfrog, copy_particle_kinetic_leapfrog
  use data_structure,     only: type_node, type_node_list, type_element, type_element_list
  use mod_parameters
  
  implicit none

  !> Data from mod_parameters for use in the k - to be removed in future using CUDA definitions in mod_parameters
  type fields_meta
     real*8               :: mu_zero, mass_proton
     real*8               :: F0, central_mass, central_density, tstep
  end type fields_meta

  !> Partial replication of fields type avoiding polymorphism for use in the kernels
  type fields_linear_device
     type(type_node_list)       :: node_list        !< Current node list
     type(type_element_list)    :: element_list     !< Current element list
     logical                    :: static           !< if true do not time interpolate
     logical                    :: flag_zero_dpsidt !< if true, P_time(1) = dpsi/dt = 0
     real*8                     :: time_now         !< Time of current restart file (SI units)
     real*8                     :: time_prev        !< Time of previous restart file (SI units)
     type(fields_meta)          :: meta
  end type fields_linear_device

  !> Partial replication of particle group type avoinding polymorphism for use in the kernels
  type particle_group_device
     integer                         :: Z
     real*8                          :: mass
     real*8                          :: dt
     type(particle_kinetic_leapfrog), managed, allocatable, dimension(:) :: particles
  end type particle_group_device
     
  
  private
  public particle_group_device, fields_linear_device, particle_kinetic_leapfrog_loop, copy_device_data

contains

  !> Principle loop to push a group of particles using CUDA kernels
  subroutine particle_kinetic_leapfrog_loop( sim  , n_steps , timestep , particle_start_time , return_particle_groups )
    type(particle_sim), intent(inout)                                    :: sim
    type(particle_group), dimension(:), allocatable                      :: return_particle_groups
    integer, intent(inout)                                               :: n_steps
    real*8, intent(inout)                                                :: particle_start_time, timestep

    type(particle_group_device), managed, dimension(:), allocatable      :: particle_groups
    type(fields_linear_device), managed, allocatable                     :: fields
    integer      :: n_groups, n_particles, tBlock_size, istat, i
    real*8       :: start_time
    type(dim3)   :: grid, tBlock

    ! Print out some info on the GPU
    if(sim%my_id==0) call device_query()

    n_steps = 1000
    timestep = 1e-10
    tBlock_size = 256
    tBlock = dim3(tBlock_size,1,1)

    ! Set the device to use for this process 
    istat = cudaSetDevice(sim%my_id)
    if (istat /= cudaSuccess) write(*,*) cudaGetErrorString(istat)

    start_time = MPI_WTime()
    call copy_device_data( sim , particle_groups , fields )
    write(*,*) "Proc ",sim%my_id," data copy completed in ",MPI_Wtime()-start_time," s"

    n_groups = size( sim%groups, 1)
    do i = 1, n_groups
       n_particles = size( sim%groups(i)%particles,1)
       grid = dim3(ceiling(real(n_particles)/tBlock_size),1,1)
       start_time = MPI_Wtime()
       write(*,*) "Group", i," launching ",n_particles, "particles"
       call particle_kinetic_leapfrog_loop_kernel<<<grid, tBlock>>>(n_particles, particle_groups(i), fields, n_steps, timestep, particle_start_time )
       istat = cudaDeviceSynchronize()
       if (istat /= cudaSuccess) write(*,*) cudaGetErrorString(istat)
       write(*,*) "Group", i," kernel completed in ",MPI_Wtime()-start_time," s"
    end do

    call copy_particle_groups_device( particle_groups, return_particle_groups )

  end subroutine particle_kinetic_leapfrog_loop


  !> The kernel to be launched to push all particles in a group
  attributes(global) subroutine particle_kinetic_leapfrog_loop_kernel( n_particles, group_particles, fields, n_steps, timestep, particle_start_time )
    use mod_boris,          only: boris_push_cylindrical
    type(particle_group_device), managed, intent(inout)         :: group_particles
    type(fields_linear_device), managed , intent(inout)         :: fields
    integer, value, intent(in)                                  :: n_particles, n_steps
    real*8, value, intent(in)                                   :: timestep, particle_start_time

    type(particle_kinetic_leapfrog)         :: particle_tmp
    integer                                 :: i,j
    real*8                                  :: t, E(3), B(3), psi, U
    real*8                                  :: rz_old(2), st_old(2)
    integer                                 :: i_elm_old, ifail

    i = threadIdx%x + (blockIdx%x-1) * blockDim%x 
    if ( i <= n_particles ) then
!       write(*,*) "+" ! simple way to count particles
       call copy_particle_kinetic_leapfrog( group_particles%particles(i) , particle_tmp )
       do j=1,n_steps
          if (particle_tmp%i_elm .le. 0) then
!             write(*,*) "-" ! simple way to count lost particles
             exit
          endif

          t = particle_start_time + (j-1)*timestep
          call calc_EBpsiU_device(fields, t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

          rz_old    = particle_tmp%x(1:2)
          st_old    = particle_tmp%st
          i_elm_old = particle_tmp%i_elm

          if (particle_tmp%i_elm .gt. 0) then
             call boris_push_cylindrical(particle_tmp, group_particles%mass, E, B, timestep)                 
             call find_rz_nearby_device(fields%node_list,fields%element_list,rz_old(1),rz_old(2), &
                  st_old(1),st_old(2),i_elm_old,particle_tmp%x(1),particle_tmp%x(2), particle_tmp%st(1), &
                  particle_tmp%st(2), particle_tmp%i_elm, ifail)                
          endif

       end do
       call copy_particle_kinetic_leapfrog( particle_tmp , group_particles%particles(i) )
    end if

  end subroutine particle_kinetic_leapfrog_loop_kernel

  !> Obtain information on the devices available, not needed for pushing
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

  !> Takes a particle sim type and copies the field and particle data into non-polymorphic
  !> types using the managed attribute that can be used in the kernels
  subroutine copy_device_data( sim , particle_groups , fields )
    type(particle_sim), intent(inout)                                               :: sim
    type(particle_group_device), managed, dimension(:), allocatable, intent(inout)  :: particle_groups
    type(fields_linear_device), managed , allocatable, intent(inout)                :: fields

    write(*,*) "Copying device data"

    call copy_particle_groups( sim%groups , particle_groups )

    call copy_fields_device( sim , fields )

  end subroutine copy_device_data

  !> Copies particle_groups to managed non-polymorphic particle_group_device groups
  subroutine copy_particle_groups( particle_groups_in, particle_groups_out )
    type(particle_group), dimension(:), intent(inout)                              :: particle_groups_in
    type(particle_group_device), managed, dimension(:), allocatable, intent(inout) :: particle_groups_out

    integer :: n_groups, n_particles, i
    
    n_groups = size(particle_groups_in)
    allocate(particle_groups_out(n_groups))

    do i=1,n_groups
       n_particles = size(particle_groups_in(i)%particles,1)
       allocate(particle_groups_out(i)%particles(n_particles))
       call copy_one_particle_group( particle_groups_in(i), particle_groups_out(i) )
    end do
    
  end subroutine copy_particle_groups

  !> Copies managed particle_group_device groups to  particle_groups
  subroutine copy_particle_groups_device( particle_groups_in, particle_groups_out )
    type(particle_group_device), managed, dimension(:), intent(inout) :: particle_groups_in
    type(particle_group), dimension(:), allocatable, intent(inout)    :: particle_groups_out

    integer :: n_groups, n_particles, i
    
    n_groups = size(particle_groups_in)
    allocate(particle_groups_out(n_groups))

    do i=1,n_groups
       n_particles = size(particle_groups_in(i)%particles,1)
       allocate(particle_kinetic_leapfrog::particle_groups_out(i)%particles(n_particles))
       call copy_one_particle_group_device( particle_groups_in(i), particle_groups_out(i) )
    end do
    
  end subroutine copy_particle_groups_device
  
  !> Copies a number of particle_groups from a particle_group to a particle_group_device type
  subroutine copy_one_particle_group(group_particles_in, group_particles_out)
    type(particle_group), intent(inout)                   :: group_particles_in
    type(particle_group_device), managed, intent(inout)   :: group_particles_out

    integer :: n_particles, i
    
    n_particles = size(group_particles_in%particles,1)

    group_particles_out%Z    = group_particles_in%Z
    group_particles_out%mass = group_particles_in%mass
    group_particles_out%dt   = group_particles_in%dt

    select type(p => group_particles_in%particles)
    type is (particle_kinetic_leapfrog)
       do i=1,n_particles
          call copy_particle_kinetic_leapfrog( p(i), group_particles_out%particles(i) )
       end do
    end select    
    
  end subroutine copy_one_particle_group
  
  !> Copies a number of particle_groups from a managed particle_group_device to a particle_group type
  subroutine copy_one_particle_group_device(group_particles_in, group_particles_out)
    type(particle_group), intent(inout)                   :: group_particles_out
    type(particle_group_device), managed, intent(inout)   :: group_particles_in

    integer :: n_particles, i
    
    n_particles = size(group_particles_in%particles,1)

    group_particles_out%Z    = group_particles_in%Z
    group_particles_out%mass = group_particles_in%mass
    group_particles_out%dt   = group_particles_in%dt

    select type(p => group_particles_out%particles)
    type is (particle_kinetic_leapfrog)
       do i=1,n_particles
          call copy_particle_kinetic_leapfrog( group_particles_in%particles(i) , p(i) )
       end do
    end select    
    
  end subroutine copy_one_particle_group_device


  subroutine copy_fields_device( sim , fields )
    type(particle_sim), intent(inout)                       :: sim
    type(fields_linear_device), managed, allocatable, intent(inout)  :: fields

    allocate(fields)

    call copy_element_list( sim%fields%element_list, fields%element_list )

    call copy_node_list( sim%fields%node_list, fields%node_list)

    call copy_fields_meta( fields%meta )

  end subroutine copy_fields_device

  subroutine copy_element_list( in , out )
    type(type_element_list), intent(in)        :: in
    type(type_element_list), intent(inout)     :: out

    integer  :: i

    out%n_elements = in%n_elements

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

  subroutine copy_node_list( in , out )
    type(type_node_list), intent(in)        :: in
    type(type_node_list), intent(inout)     :: out

    integer  :: i

    out%n_nodes = in%n_nodes
    out%n_dof = in%n_dof

    do i=1,n_nodes_max
       call copy_node( in%node(i) , out%node(i) )
    end do

  end subroutine copy_node_list

  !> Copy a node type from the host to a GPU device
  subroutine copy_node( in , out ) 

    type(type_node), intent(in)  :: in 
    type(type_node), intent(out) :: out 

    out%x              = in%x
    out%values         = in%values
    out%deltas         = in%deltas
#ifdef fullmhd
    out%psi_eq         = in%psi_eq
    out%Fprof_eq       = in%Fprof_eq
#elif altcs
    out%psi_eq         = in%psi_eq
#endif
    out%index          = in%index
    out%boundary       = in%boundary
    out%boundary_index = in%boundary_index
    out%axis_node      = in%axis_node
    out%parents        = in%parents
    out%parent_elem    = in%parent_elem
    out%ref_lambda     = in%ref_lambda
    out%ref_mu         = in%ref_mu
    out%constrained    = in%constrained

  end subroutine copy_node


  !> Copy meta data from the element_list and node_list types into a fields_meta type
  !> for use in the interpolation routines
  subroutine copy_fields_meta( meta )

    use constants, only: mu_zero, mass_proton
    use phys_module, only: F0, tstep, central_mass, central_density

    type(fields_meta), intent(inout) :: meta

    meta%mu_zero          = mu_zero
    meta%mass_proton      = mass_proton     
    meta%F0               = F0     
    meta%central_mass     = central_mass     
    meta%central_density  = central_density    
    meta%tstep            = tstep

  end subroutine copy_fields_meta



  !> Version of routine from mod_fields for running on a GPU device
  attributes(device) subroutine calc_EBpsiU_device(fields, time, i_elm, st, phi, E, B, psi, U)

    type(fields_linear_device), managed, intent(in)    :: fields
    real*8, intent(in)  :: time
    integer, intent(in) :: i_elm !< JOREK element index
    real*8, intent(in)  :: st(2) !< element-local coordinates
    real*8, intent(in)  :: phi !< toroidal angle
    real*8, intent(out) :: E(3) !< Electric field [V/m]
    real*8, intent(out) :: B(3) !< Magnetic field [T]
    real*8, intent(out) :: psi !< psi in JOREK units
    real*8, intent(out) :: u !< velocity stream function in m/s
    ! Internal parameters
    !    integer, parameter :: i_var(2) = [1,2]
    integer             :: i_var(2) 
    real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_time(2) ! Placeholder for evaluating variables and derivatives locally
    ! Values
    real*8             :: R, R_s, R_t, Z, Z_s, Z_t
    ! Others
    real*8             :: inv_st_jac, R_inv
    real*8             :: psi_R, psi_Z, U_R, U_Z, U_phi, t_norm

    i_var(1) = 1
    i_var(2) = 2

    t_norm  = sqrt(fields%meta%mu_zero * fields%meta%mass_proton * fields%meta%central_mass * fields%meta%central_density * 1.d20) ! 1 jorek time unit in seconds

    ! Interpolate the fields to get psi and U at the current position (and the
    ! changes u_n - u(n-1))
    call do_interp_PRZ_device(fields, time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

    R_inv = 1.d0/R
    inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)

    ! Calculate the derivatives to R and Z
    psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
    psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
    U_R      = (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
    U_Z      = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
    U_phi    = P_phi(2)

    ! Update psi and U
    psi = P(1)
    U   = P(2)/t_norm

    ! Set dpsi/dt to 0 if flag is true
    if(fields%flag_zero_dpsidt) P_time(1) = 0.d0

    ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
    B     = [+psi_Z, -psi_R, fields%meta%F0] * R_inv

    ! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
    ! See http://jorek.eu/wiki/doku.php?id=u_phi
    E     = [-fields%meta%F0*U_R, -fields%meta%F0*U_Z, -fields%meta%F0*U_phi*R_inv]/t_norm
    E(3)  = E(3) - R_inv*P_time(1) ! because this is not normalized with t_norm

  end subroutine calc_EBpsiU_device

  !> Interpolate a variable at a specific position (with phi), with first derivatives only
  attributes(device) subroutine do_interp_PRZ_device(fields, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
    use mod_linear 
    
    type(fields_linear_device), managed, intent(in)  :: fields
    real*8,                   intent(in)           :: time !< Time at which to calculate this variable
    integer,                  intent(in)           :: i_elm
    integer,                  intent(in)           :: n_v, i_v(n_v)
!!$    integer,                  intent(in)           :: n_v, i_v(2)
    real*8,                   intent(in)           :: s, t, phi
    real*8,                   intent(out)          :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v), P_time(n_v)
!!$    real*8,                   intent(out)          :: P(2), P_s(2), P_t(2), P_phi(2), P_time(2)
    real*8,                   intent(out)          :: R, R_s, R_t, Z, Z_s, Z_t

    real*8                 :: df, dt
!!$  real*8, dimension(n_v) :: Pd, Pd_s, Pd_t, Pd_phi
    real*8, dimension(2) :: Pd, Pd_s, Pd_t, Pd_phi
    real*8                 :: t_jorek

    ! JOREK time step in seconds
    t_jorek = fields%meta%tstep*sqrt(fields%meta%mu_zero * fields%meta%mass_proton * &
         fields%meta%central_mass * fields%meta%central_density * 1.d20)
    P_time = 0.d0

    !> interpolate values
    call interp_PRZ_device(fields%node_list, fields%element_list,i_elm ,i_v ,n_v ,s ,t , phi, &
         P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t, .False.)

    !> interpolate differentials
    if(t_jorek .gt. 0.d0) then
       call interp_PRZ_device(fields%node_list, fields%element_list, i_elm, i_v, n_v, s, t, phi, &
            Pd, Pd_s, Pd_t, Pd_phi, R, R_s, R_t, Z, Z_s, Z_t, .True.)
       if(abs(fields%time_now-fields%time_prev) .gt. 1d-10 .and. .not. fields%static) then
          !> compute time fraction df
          dt = 1.d0/(fields%time_now - fields%time_prev)
          df = (fields%time_now - time)*dt
          !> apply linear interpolation
          P     = linear_interp_differentials(n_v,P,Pd,df)
          P_s   = linear_interp_differentials(n_v,P_s,Pd_s,df)
          P_t   = linear_interp_differentials(n_v,P_t,Pd_t,df)
          P_phi = linear_interp_differentials(n_v,P_phi,Pd_phi,df)
       else
          dt = 1.d0/t_jorek
       endif
       !> compute time derivative
       P_time = linear_interp_differentials_dt(n_v,Pd,dt) 
    endif

  end subroutine do_interp_PRZ_device



  !> This subroutine interpolates some variables at a specific position within one element at a given position (s,t)
  attributes(device) subroutine interp_PRZ_device( node_list, element_list, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t, deltas)
    use mod_basisfunctions
    use mod_interp,         only: sincosperiod_moivre

    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list    
    integer, intent(in)                            :: i_elm
    integer, intent(in)                            :: n_v, i_v(2)
!!$    integer, intent(in)                            :: n_v, i_v(n_v)
    real*8, intent(in)                             :: s, t, phi
!!$    real*8, intent(out)                            :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v)
    real*8, intent(out)                            :: P(2), P_s(2), P_t(2), P_phi(2)
    real*8, intent(out)                            :: R, R_s, R_t, Z, Z_s, Z_t
    logical, intent(in)                            :: deltas

    ! --- Local variables
    real*8  :: H(4,4), H_s(4,4), H_t(4,4), HZ(n_tor), dHZ(n_tor)
    integer :: kv, iv, kf, i
    real*8  :: values(n_tor,n_order+1,2,n_vertex_max) ! replaced n_v with 2
!!$    real*8  :: values(n_tor,n_order+1,n_v,n_vertex_max)
    real*8  :: xR(n_order+1,n_vertex_max), xZ(n_order+1,n_vertex_max)
    real*8  :: sizes(n_order+1), v, vp

    ! 7% exec time
    call basisfunctions_T(s,t,H,H_s,H_t)

    P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0

    ! 7% exec time
    call sincosperiod_moivre(phi, HZ, dHZ)

    ! 30% exec time
    ! Preload values and premultiply with sizes(:,kv)
    do kv = 1,n_vertex_max
       iv = element_list%element(i_elm)%vertex(kv)
       sizes(:) = element_list%element(i_elm)%size(kv,:)
       if (deltas) then
          do i = 1, n_v
             do kf=1,n_order+1
                values(1:n_tor,kf,i,kv) = node_list%node(iv)%deltas(1:n_tor,kf,i_v(i)) * sizes(kf)
             end do
          end do
       else
          do i = 1, n_v
             do kf=1,n_order+1
                values(1:n_tor,kf,i,kv) = node_list%node(iv)%values(1:n_tor,kf,i_v(i)) * sizes(kf)
             end do
          end do
          do i = 1, n_v
             do kf=1,4
                values(1,kf,i,kv) = node_list%node(iv)%values(1,kf,i_v(i)) * sizes(kf)
             end do
          end do
       end if
       xR(:,kv) = node_list%node(iv)%x(1,:,1) * sizes(:)
       xZ(:,kv) = node_list%node(iv)%x(1,:,2) * sizes(:)
    enddo

    ! together 7%
    R   = sum(xR*H)
    R_s = sum(xR*H_s)
    R_t = sum(xR*H_t)
    Z   = sum(xZ*H)
    Z_s = sum(xZ*H_s)
    Z_t = sum(xZ*H_t)

    ! 40% exec time
    do kv = 1, n_vertex_max
       do i = 1, n_v
          do kf = 1, n_order+1
             v = dot_product(values(1:n_tor,kf,i,kv),HZ(1:n_tor))
             P(i)     = P(i)     + v * H(kf, kv)
             P_s(i)   = P_s(i)   + v * H_s(kf, kv)
             P_t(i)   = P_t(i)   + v * H_t(kf, kv)
             vp = dot_product(values(1:n_tor,kf,i,kv),dHZ(1:n_tor))
             P_phi(i) = P_phi(i) + vp * H(kf, kv)
          enddo
       enddo
    enddo

  end subroutine interp_PRZ_device

  attributes(device) subroutine find_RZ_nearby_device(node_list, element_list, R_old, Z_old, s_old, t_old, i_elm_old, &
       R_new, Z_new, s_new, t_new, i_elm_new, ifail)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use mod_neighbours, only : coord_in_neighbour
    use mod_find_rz_nearby, only : try_interp
    implicit none
    !> Input parameters
    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8,                   intent(in)    :: R_old, Z_old !< The old R,Z location
    real*8,                   intent(in)    :: R_new, Z_new !< The new R,Z location
    real*8,                   intent(in)    :: s_old, t_old !< The old st location (used to compute a guess)
    integer,                  intent(in)    :: i_elm_old
    real*8,                   intent(out)   :: s_new, t_new !< The found new coordinates
    integer,                  intent(out)   :: i_elm_new
    integer,                  intent(out)   :: ifail !< if ifail = -1 the position could not be found in the grid.
    !< ifail > 0 indicates various other cases

    !> Accuracy defaults (tolerances are squared!, units of element size)
    real*8,  parameter :: element_tolerance   = 1.d-24 !< Tolerance for finding a position inside an element
    integer, parameter :: newton_iter_max     = 8 !< Number of iterations to try

    !> Internal variables
    integer :: newton_iter_number, i_elm_tmp
    real*8 :: inv_st_jac_det, R_s, R_t, Z_s, Z_t
    real*8 :: st_step(2), x_step(2), x_tmp(2), st_new(2), x_new(2) ! x_step = (R,Z) of trial position
    real*8 :: err2, err2_old, dist(2), fact

    ! Check if element is valid
    if (i_elm_old .lt. 1 .or. i_elm_old .gt. element_list%n_elements) then
       call find_RZ_device(node_list,element_list,R_new,Z_new,x_step(1),x_step(2),i_elm_new,s_new,t_new,ifail)
       return
    end if
    ! Setup initial values
    x_step = [R_old,Z_old] ! start at the current position
    i_elm_new = i_elm_old ! start in the current element
    st_new = [s_old,t_old] ! start at the old position
    x_new = [R_new,Z_new]
    ! Find the jacobian at the current s and t position
    call try_interp(node_list,element_list,i_elm_new,st_new,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
    err2 = dot_product(x_step-x_new,x_step-x_new)
    ifail=0


    ! Newton iteration to find s and t in or out of this element
    do newton_iter_number = 1, newton_iter_max
       ! Perform newton iteration by calculating the inverse of the jacobian matrix explicitly
       err2_old = err2

       ! Calculate the trial newton step
       st_step(1) = ( Z_t * (x_new(1)-x_step(1)) - R_t * (x_new(2)-x_step(2))) * inv_st_jac_det
       st_step(2) = (-Z_s * (x_new(1)-x_step(1)) + R_s * (x_new(2)-x_step(2))) * inv_st_jac_det

       ! Limit this step if it goes outside of the element
       dist = merge(1-st_new,st_new,st_step .gt. 0) ! dist = 1-s if step > 0, s if step < 0 (distance to 0 or 1)
       fact = maxval(abs(st_step)/max(dist,1d-30)) ! if fact>=1 we are on the boundary
       ! (it is the overshoot: i.e. how many times we overshoot the boundary with one st_step)

       if (fact .ge. 1.d0-1d-12) then
          st_new = st_new + st_step/fact
#ifdef DEBUG
          call try_interp(node_list,element_list,i_elm_new,st_new,x_tmp,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
#endif
          i_elm_tmp = i_elm_new
          call coord_in_neighbour(node_list,element_list,i_elm_tmp,i_elm_new,st_new)
          if (i_elm_new .lt. 0) then
             call find_RZ_device(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,s_new,t_new,ifail)
             if (ifail .ne. 0) i_elm_new = 0
          end if
          if (i_elm_new .eq. 0) then ! No element on that side, particle is lost
             i_elm_new = - i_elm_tmp ! Save position of particle
             ! Calculate new R and Z in x_new
             call try_interp(node_list,element_list,i_elm_tmp,st_new,x_new,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
             ! Set new element-local coordinates for the point on the axis
             s_new = st_new(1)
             t_new = st_new(2)
             ifail = -1
             return
          end if

          call try_interp(node_list,element_list,i_elm_new,st_new,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
#ifdef DEBUG 
          if (norm2(x_step-x_tmp) .gt. 1d-8) then
             !write(*,*) "ERROR on element edge crossing", x_step, x_tmp, norm2(x_step-x_tmp), &
             !i_elm_new, i_elm_old, i_elm_tmp, "deleting particle"
             !call exit(1)
             i_elm_new = 0
             return
          end if
#endif
       else
          st_new = st_new + st_step
          call try_interp(node_list,element_list,i_elm_new,st_new,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
       end if
       err2 = dot_product(x_step-x_new,x_step-x_new)
       s_new = st_new(1)
       t_new = st_new(2)

       if (err2 < element_tolerance) exit
    enddo


    if (ieee_is_nan(err2)) then
       !write(*,*) "WARNING: NaN encountered after newton iteration, using find_RZ"
       call find_RZ_device(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,s_new,t_new,ifail)
       if (ifail .eq. 0) ifail=2
       return
    endif
    if (newton_iter_number .gt. newton_iter_max) then
       !write(*,"(A,i4,A,i5,A,2g14.6,A,3g14.6)") "WARNING: iteration for st did not converge after", newton_iter_max, " tries in element ", i_elm_new, &
       !" using find_RZ", x_new, "err2(old)/convergence: ", err2, err2_old, err2_old/err2
       !write(*,"(A,2g16.8)") "Find_RZ at ", x_new
       call find_RZ_device(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,s_new,t_new,ifail)
       if (ifail .eq. 0) ifail=3
       return
    endif
  end subroutine find_RZ_nearby_device


  !> Auxiliary subroutine for find_RZ_nearby
  attributes(device) pure subroutine try_interp_device(node_list,element_list,i_elm,st,x,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
    use mod_interp, only : interp_RZ

    !> Input parameters
    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    real*8,                   intent(in)    :: st(2)
    integer,                  intent(in)    :: i_elm
    real*8,                   intent(out)   :: x(2), R_s, R_t, Z_s, Z_t, inv_st_jac_det
    real*8 :: jac

    call interp_RZ(node_list,element_list,i_elm,st(1),st(2),x(1),R_s,R_t,x(2),Z_s,Z_t)
    ! Guard against the determinant being close to zero
    jac = R_s * Z_t - R_t * Z_s
    if (abs(jac) .lt. 1d-8) then
       inv_st_jac_det = sign(1d8, jac) 
    else
       inv_st_jac_det = 1.d0/(jac)
    end if
  end subroutine try_interp_device

  attributes(device) subroutine find_RZ_device(node_list,element_list,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
    !-------------------------------------------------------------------------
    !< Find all elements for which minmax is correct and run find_RZ_single on those.
    !< Return the first result.
    !-------------------------------------------------------------------------
    implicit none

    type (type_node_list), intent(in)    :: node_list
    type (type_element_list), intent(in) :: element_list
    real*8, intent(in)     :: R_find, Z_find
    real*8, intent(out)    :: R_out,Z_out,s_out,t_out
    integer, intent(inout) :: ielm_out
    integer, intent(out)   :: ifail

    integer :: k, ielm_in
    !integer, dimension(:), allocatable :: i_elms
    integer, dimension(20) :: i_elms

    ielm_in = ielm_out
    ielm_out = 0
    !call elements_containing_point(R_find, Z_find, i_elms)
    !HJL Brute force to avoid rtree code
    !allocate(i_elms(4))
    do k=1, 4
       i_elms(4*k+1:4*k+4) = element_list%element(i_elms(k))%neighbours
    enddo


    ! then loop through all
    do k=1,20 !size(i_elms)
       call find_RZ_single_device(node_list,element_list,i_elms(k),R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
       if (ifail .eq. 0) exit
    enddo

    if (ielm_out .eq. 0) ifail = 99
    if (ifail .eq. 999) ielm_out = 0 ! Otherwise testing ielm=0 on output does not
    ! work anymore (and we don't always check ifail)

  end subroutine find_RZ_device

  attributes(device) subroutine find_RZ_single_device(node_list,element_list,i_elm,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
    !-------------------------------------------------------------------------
    !< solves two non-linear equations using Newtons method (from numerical recipes)
    !< LU decomposition replaced by explicit solution of 2x2 matrix.
    !<
    !< finds the crossing of two coordinate lines given as a series of cubics in element
    !< i_elm
    !-------------------------------------------------------------------------
    use mod_interp, only : interp_RZ
    implicit none

    type (type_node_list), intent(in)    :: node_list
    type (type_element_list), intent(in) :: element_list
    integer, intent(in)    :: i_elm
    real*8, intent(in)     :: R_find, Z_find
    real*8, intent(out)    :: R_out,Z_out,s_out,t_out
    integer, intent(out)   :: ielm_out
    integer, intent(out)   :: ifail

    integer :: i, ntrial, istart
    real*8  :: RRg1,dRRg1_dr,dRRg1_ds
    real*8  :: ZZg1,dZZg1_dr,dZZg1_ds
    real*8  :: tolx, tolf, errx, errf, temp, dis
    real*8  :: x(2), FVEC(2), FJAC(2,2), p(2)

    ntrial = 20
    tolx = 1.d-8
    tolf = 1.d-15

    ielm_out = i_elm ! Since we only test a single element

    do istart = 1,5

       if (istart .eq. 1) then
          x(1) = 0.5d0
          x(2) = 0.5d0
       elseif (istart .eq. 2) then
          x(1) = 0.75d0
          x(2) = 0.75d0
       elseif (istart .eq. 3) then
          x(1) = 0.75d0
          x(2) = 0.25d0
       elseif (istart .eq. 4) then
          x(1) = 0.25d0
          x(2) = 0.75d0
       elseif (istart .eq. 5) then
          x(1) = 0.25d0
          x(2) = 0.25d0
       endif

       ifail = 999

       do i=1,ntrial
          !HJL
              call interp_RZ(node_list,element_list,i_elm,x(1),x(2),RRg1,dRRg1_dr,dRRg1_ds, &
                                                              ZZg1,dZZg1_dr,dZZg1_ds)
          FVEC(1)   = RRg1 - R_find
          FVEC(2)   = ZZg1 - Z_find
          FJAC(1,1) = dRRg1_dr
          FJAC(1,2) = dRRg1_ds
          FJAC(2,1) = dZZg1_dr
          FJAC(2,2) = dZZg1_ds

          errf=abs(fvec(1))+abs(fvec(2))

          !      write(*,'(A,i3,8e16.8)') ' newton   : ',i,errf,errx,x,RRg1,R_find,ZZg1,Z_find
          !      write(*,'(A,i3,8e16.8)') ' newton   : ',i,dRRg1_dr,dRRg1_ds,dZZg1_dr,dZZg1_ds

          if (errf .le. tolf) then

             s_out     = x(1)
             t_out     = x(2)

             ielm_out  = i_elm
             R_out     = RRg1
             Z_out     = ZZg1

             !        write(*,'(A,i3,4e16.8)') ' newton (1) : ',i,errf,errx,x

             ifail = 0
             return
          endif

          p = -fvec

          temp = p(1)
          dis  = fjac(2,2)*fjac(1,1)-fjac(1,2)*fjac(2,1)

          if (dis .ne. 0.d0) then
             p(1) = (fjac(2,2)*p(1)-fjac(1,2)*p(2))/dis
             p(2) = (fjac(1,1)*p(2)-fjac(2,1)*temp)/dis
          else
             exit
          endif

          errx=abs(p(1)) + abs(p(2))

          p = min(p,+0.25d0)
          p = max(p,-0.25d0)

          x = x + p

          x = max(x,+0.d0)
          x = min(x,+1.d0)

          if (errx .le. tolx) then

             s_out     = x(1)
             t_out     = x(2)

             ielm_out  = i_elm
             R_out     = RRg1
             Z_out     = ZZg1

             !        write(*,'(A,i3,4e16.8)') ' newton (2) : ',i,errf,errx,x

             ifail = 0
             return
          endif

       enddo
    enddo
  end subroutine find_RZ_single_device

 end module mod_particle_kernels
#endif
