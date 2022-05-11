#ifdef CUDA_KERNELS
!> Contains routines for offloading particle computation to GPU devices using CUDA kernels
!> At present only the particle_kinetic_leapfrog type is supported

module mod_particle_kernels

  use cudafor
  use mpi
  use mod_particle_sim,   only: particle_sim
  use mod_particle_types, only: particle_kinetic_leapfrog, copy_particle_kinetic_leapfrog
  use mod_boris,          only: boris_push_cylindrical
  use mod_basisfunctions
  use mod_interp
  use data_structure
  use mod_linear
  
  implicit none

  type fields_meta
     real*8               :: mu_zero, mass_proton
     real*8               :: F0, central_mass, central_density, tstep
!!$     integer              :: n_tor, n_period, n_vertex_max, n_order
!!$     real*8               :: time_now, time_prev
  end type fields_meta
  
  type fields_linear_device
     type(type_node_list)       :: node_list        !< Current node list
     type(type_element_list)    :: element_list     !< Current element list
     logical                    :: static           !< if true do not time interpolate
     logical                    :: flag_zero_dpsidt !< if true, P_time(1) = dpsi/dt = 0
     real*8                     :: time_now         !< Time of current restart file (SI units)
     real*8                     :: time_prev        !< Time of previous restart file (SI units)
     type(fields_meta)          :: meta
  end type fields_linear_device

  private
  public fields_linear_device, particle_kinetic_leapfrog_loop, copy_device_data
  
contains


  subroutine particle_kinetic_leapfrog_loop( sim  , n_steps , timestep , particle_start_time )
    type(particle_sim), intent(inout)                                    :: sim
    integer, intent(inout)                                               :: n_steps
    real*8, intent(inout)                                                :: particle_start_time, timestep
    
    type(particle_kinetic_leapfrog), managed, dimension(:), allocatable  :: particles
    type(fields_linear_device), managed, allocatable                     :: fields
    integer      :: n_particles, tBlock_size, istat
    real*8       :: start_time
    type(dim3)   :: grid, tBlock

    ! Print out some info on the GPU
    if(sim%my_id==0) call device_query()

    n_steps = 1000
    timestep = 1e-10
    n_particles = size( sim%groups(1)%particles,1)
    tBlock_size = 256
    tBlock = dim3(tBlock_size,1,1)
    grid = dim3(n_particles/tBlock_size,1,1)

    istat = cudaSetDevice(sim%my_id)
    if (istat /= cudaSuccess) write(*,*) cudaGetErrorString(istat)

    start_time = MPI_WTime()
    call copy_device_data( sim , particles , fields )
    write(*,*) "Proc ",sim%my_id," data copy completed in ",MPI_Wtime()-start_time," s"

    start_time = MPI_WTime()
    call particle_kinetic_leapfrog_loop_kernel<<<grid, tBlock>>>(n_particles, particles, fields, n_steps, timestep, particle_start_time )
    istat = cudaDeviceSynchronize()
    if (istat /= cudaSuccess) write(*,*) cudaGetErrorString(istat)
    write(*,*) "Proc ",sim%my_id," kernel completed in ",MPI_Wtime()-start_time," s"
   
  end subroutine particle_kinetic_leapfrog_loop

  attributes(global) subroutine particle_kinetic_leapfrog_loop_kernel( n_particles, particles, fields, n_steps, timestep, particle_start_time )
    type(particle_kinetic_leapfrog), managed, dimension(:), intent(inout)  :: particles
    type(fields_linear_device), managed , intent(inout)                    :: fields
    integer, value, intent(in)                                             :: n_particles, n_steps
    real*8, value, intent(in)                                              :: timestep, particle_start_time

    type(particle_kinetic_leapfrog)         :: particle_tmp
    integer                                 :: i,j
    real*8                                  :: t, E(3), B(3), psi, U, mass
    real*8                                  :: rz_old(2), st_old(2)
    integer                                 :: i_elm_old, ifail

    mass = 2
    i = threadIdx%x + (blockIdx%x-1) * blockDim%x 
    if ( i <= n_particles ) then
       call copy_particle_kinetic_leapfrog( particles(i) , particle_tmp )
       do j=1,n_steps
          if (particle_tmp%i_elm .le. 0) then
             write(*,*) "Losing particle",j
             exit
          endif

          t = particle_start_time + (j-1)*timestep
          call calc_EBpsiU_device(fields, t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

          rz_old    = particle_tmp%x(1:2)
          st_old    = particle_tmp%st
          i_elm_old = particle_tmp%i_elm

          if (particle_tmp%i_elm .gt. 0) then
             call boris_push_cylindrical(particle_tmp, mass, E, B, timestep)                 
          endif

       end do
       call copy_particle_kinetic_leapfrog( particle_tmp , particles(i) )
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
    type(fields_linear_device), managed , allocatable, intent(inout)                   :: fields
    
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
    
!    type(particle_sim), intent(in) :: sim
    type(fields_meta), intent(inout) :: meta

!!$    fields_meta%static           = sim%fields%static
!!$    fields_meta%flag_zero_dpsidt = sim%fields%flag_zero_dpsidt
!!$    fields_meta%n_elements       = sim%fields%element_list%n_elements
!!$    fields_meta%n_nodes          = sim%fields%node_list%n_nodes
!!$    fields_meta%n_dof            = sim%fields%node_list%n_dof

! Jorek parameters
!!$    meta%n_period         = n_period
!!$    meta%n_tor            = n_tor
!!$    meta%n_vertex_max     = n_vertex_max
!!$    meta%n_order          = n_order

    meta%mu_zero          = mu_zero
    meta%mass_proton      = mass_proton     
    meta%F0               = F0     
    meta%central_mass     = central_mass     
    meta%central_density  = central_density    
    meta%tstep            = tstep
    
  end subroutine copy_fields_meta


  
  !> Version of routine from mod_fields for running on a GPU device
  attributes(device) subroutine calc_EBpsiU_device(fields, time, i_elm, st, phi, E, B, psi, U)

    type(fields_linear_device), intent(in)    :: fields
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

    type(fields_linear_device), intent(in)         :: fields
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

end module mod_particle_kernels
#endif
