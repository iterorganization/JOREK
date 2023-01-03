!> Read a set of particles and static fields and calculate the fieldline trace
!> Command-line arguments: part_restart.h5 < inputfile
!>
!> TODO: * Turn this into a diagnostic type and enable as a default step for simulations
program calc_phi_zero
use constants
use mod_particle_sim
use mod_particle_types
use mod_io_actions
use mod_fields
use mod_fields_linear
use mpi
use mod_event
use mod_boris
use mod_fieldline_euler
use hdf5_io_module
use mod_find_rz_nearby
implicit none

type(particle_sim) :: sim
real*8, parameter :: dt = 1d-8
real*8, parameter :: v = 1d6
type(event) :: fieldreader, partreader
integer :: i, j, k, i_elm_old, ifail, ierr, i_elm, dir
integer(HID_T) :: file_id, group_id
logical :: link_exists
real*8 :: psi_axis, R_axis, Z_axis, s_axis, t_axis
real*8 :: E(3), B(3), rz_old(2), st_old(2), phi, phi_old, theta, theta_old, psi, U
real*8, allocatable, dimension(:,:) :: phi_zero, phi_zero_dist
type(particle_fieldline), allocatable, dimension(:) :: particles
character(len=500) :: part_file
character(len=20) :: group_name
part_file = ''

! Start up MPI, jorek
call sim%initialize(num_groups=0)
call get_command_argument(1, part_file)

partreader = event(read_action(filename=trim(part_file)))
call with(sim, partreader) 
fieldreader = event(read_jorek_fields_interp_linear(i=-1))
call with(sim, fieldreader)
! For some reason the code will segfault if you use the same event for field and particle reading
! It is also important to read the fields second. Something tricky is going on with the particle reader
! which ruins the fields that have been read because on its own this works fine.
! This is a very bad sign but I have no time to fix it now. Leave the ordering above as-is.

! Find the axis to use for calculating theta (particle is on outer midplane if theta \approx 0 and R >= R_axis)
call find_axis(0, sim%fields%node_list, sim%fields%element_list, psi_axis, R_axis, Z_axis, i_elm, s_axis, t_axis, ifail)


! Open an (existing) HDF5 file to save this data. Use a format similar to the particle output file (/groups/001/phi_zero, /groups/001/L)
call HDF5_open_or_create(part_file, file_id=file_id, ierr=ierr)
! Create group to write particle groups in if it does not exist yet
call h5lexists_f(file_id, "/groups", link_exists, ierr)
if (.not. link_exists) then
  call h5gcreate_f(file_id, "/groups", group_id, ierr)
  call h5gclose_f(group_id, ierr)
end if



! Trace each fieldline particle in both directions until we get to the outboard midplane (theta = 0)
! and take the one which has the shortest distance
do i=1,size(sim%groups,1)
  write(group_name,'(A,i0.3,A)') '/groups/', i, '/'
  allocate(     phi_zero(size(sim%groups(i)%particles),2))
  allocate(phi_zero_dist(size(sim%groups(i)%particles),2))
  allocate(particles(size(sim%groups(i)%particles)))
  do dir=1,2
    select type (p => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      do j=1,size(p)
        if (p(j)%i_elm .le. 0) then
          particles(j) = p(j)
        else
          ! Calculate magnetic field to get GC coordinate
          call sim%fields%calc_EBpsiU(0.d0, p(j)%i_elm, &
              p(j)%st, p(j)%x(3), E, B, psi, U)
          particles(j) = kinetic_leapfrog_to_gc(sim%fields%node_list, sim%fields%element_list, p(j), E, B, sim%groups(i)%mass, dt=0.d0)
          ! dt above is not per-se the right dt for this particle (since we read it from a file). Use 0 instead
        end if
      end do
    type is (particle_gc)
      do j=1,size(p)
        particles(j) = p(j)
      end do
    end select ! Add types here as needed

    do j=1,size(particles)
      particles(j)%v = v*(real(dir*2-3))
    end do

    phi_zero = 0.d0
    phi_zero_dist = 1d99 ! so we don't select it by default
    !$omp parallel do default(private) &
    !$omp shared(sim, dir, i, phi_zero, phi_zero_dist, R_axis, Z_axis, particles)
    do j=1,size(particles,1)
      if (particles(j)%i_elm .le. 0) cycle
      ! Do a single euler step forward to setup the adams-bashforth method
      call sim%fields%calc_EBpsiU(0.d0, particles(j)%i_elm, &
        particles(j)%st, particles(j)%x(3), E, B, psi, U)
      rz_old    = particles(j)%x(1:2)
      st_old    = particles(j)%st
      i_elm_old = particles(j)%i_elm
      call fieldline_euler_push_cylindrical(particles(j), B, dt)
      call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
          particles(j)%x(1), particles(j)%x(2), particles(j)%st(1), particles(j)%st(2), particles(j)%i_elm, ifail)
      if (particles(j)%i_elm .le. 0) cycle
      particles(j)%B_hat_prev = B/norm2(B)

      do k=1,200*nint(1.d0/(v*dt)) ! maximum number of steps from maximum length = q*circumference/v/dt \approx 100/vdt
        call sim%fields%calc_EBpsiU(0.d0, particles(j)%i_elm, &
          particles(j)%st, particles(j)%x(3), E, B, psi, U)
        rz_old    = particles(j)%x(1:2)
        phi_old   = particles(j)%x(3)
        st_old    = particles(j)%st
        i_elm_old = particles(j)%i_elm
        theta_old = atan2(particles(j)%x(2)-Z_axis, particles(j)%x(1)-R_axis)

        call fieldline_adams_bashforth_push_cylindrical(particles(j), B, dt)
        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
            particles(j)%x(1), particles(j)%x(2), particles(j)%st(1), particles(j)%st(2), particles(j)%i_elm, ifail)
        if (particles(j)%i_elm .le. 0) exit

        ! Check if the new position is at the outer midplane
        theta = atan2(particles(j)%x(2)-Z_axis, particles(j)%x(1)-R_axis)
        ! This is true if theta has a different sign from theta_old and they are both smaller than pi/2
        if (sign(1.d0,theta) .ne. sign(1.d0,theta_old) .and. abs(theta) .lt. 0.5d0*PI .and. abs(theta_old) .lt. 0.5d0*PI) then 
          ! interpolate linearly (small-angle approx) to find a more precise value of phi
          phi = particles(j)%x(3)
          phi_zero(j,dir) = phi - theta*(phi - phi_old)/(theta - theta_old)
          ! Calculate the distance more precisely as well
          phi_zero_dist(j,dir) = (real(k-1) + abs(theta)/abs(theta-theta_old))*v*dt
          exit
        end if
      end do ! steps
    end do ! particles
    !$omp end parallel do
  end do ! dir
  ! Select the shortest distance to the outer midplane
  ! Do this by copying into phi_zero(:,1) values from phi_zero(:,2) where phi_zero_dist(:,1) > phi_zero_dist(:,2)
  where (phi_zero_dist(:,1) .gt. phi_zero_dist(:,2)) phi_zero(:,1) = phi_zero(:,2)
  where (phi_zero_dist(:,1) .gt. phi_zero_dist(:,2)) phi_zero_dist(:,1) = -phi_zero_dist(:,2)
  call HDF5_array1D_saving(file_id, phi_zero(:,1), size(phi_zero,1),trim(group_name)//"Phi_zero")
  call HDF5_array1D_saving(file_id, phi_zero_dist(:,1), size(phi_zero_dist,1),trim(group_name)//"L")
  deallocate(phi_zero,phi_zero_dist,particles)
  write(*,*) "Written group ", i
end do ! groups
call h5fclose_f(file_id, ierr)
call h5close_f(ierr)
call sim%finalize
end program calc_phi_zero
