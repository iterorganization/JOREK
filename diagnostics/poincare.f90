!< Program for creating Poincare plots using the particle tracer.
!<
!< Usage:
!< ./poincare < input_namelist
!<
!< This tool uses either the field line or relativistic guiding center (RK4 in both cases) pushers to create Poincare plots.
!< Simulation settings are given in a file named "pncr" and it should have the following format:
!<
!< # N_tracers
!< 100
!< # n_turns
!< 200
!< # tstep in [s] for particles and [m] for field lines
!< 1.0e-8
!< # mass [AMU] (use 0.0 to trace field lines)
!< 0.00054
!< # charge [e]
!< -1
!< # pitch
!< 0.99
!< # energy [eV]
!< 1.0e3
!<
!< The Poincare data is collected using jorek_restart.h5 as the field input.
!<
program poincare

use particle_tracer
use mod_particle_io
use mod_event
use mod_fields_linear
use mod_fieldline_euler
use mod_gc_relativistic
use hdf5_io_module
use phys_module, only: xcase, xpoint
use equil_info
use mod_boundary, only: boundary_from_grid
use nodes_elements
use mod_interp, only: interp
use mod_coordinate_transforms, only: cylindrical_to_cartesian

implicit none

character(len=40) :: fnout !< File where the output is written

real*8, dimension(:,:), allocatable    :: rvals, zvals, phivals, psivals !< Coordinate arrays (nprt,norb).
integer*4, dimension(:,:), allocatable :: pncrid   !< Flag indicating which plane was crossed (nprt,norb).
integer*4, dimension(:), allocatable   :: nextslot !< Next free slot in coordinate arrays.
integer*4, dimension(:), allocatable   :: ncross   !< Number of poloidal turns
real*8, dimension(:), allocatable      :: mileage  !< How far marker has travelled in secons (in meters for field line)
real*8, dimension(:,:), allocatable    :: milvals  !< mileage but recorded separately at each crossing.

type(event)       :: fieldreader
class(*), pointer :: p
real*8    :: tstep
integer*4 :: norb, ndata, nprt, finished
integer*4 :: ifail, iprt
real*8    :: zprev, phiprev, raxis, zaxis, xyz(3), xyz0(3)
real*8    :: taxis, saxis, psiaxis, ielmaxis
real*8    :: mass, pitch, energy
integer*4 :: charge
character(len=512) :: s

! For CPU time
real*8 :: t0, t1

! Simulation options
open(21, file='pncr', status='old', action='read', iostat=ifail)

if ( ifail == 0 ) then ! pncr file exists, use it.
   read(21, '(a)') s ! read comment line (ignored)
   read(21,*) nprt
   read(21, '(a)') s
   read(21,*) norb
   read(21, '(a)') s
   read(21,*) tstep
   read(21, '(a)') s
   read(21,*) mass
   read(21, '(a)') s
   read(21,*) charge
   read(21, '(a)') s
   read(21,*) pitch
   read(21, '(a)') s
   read(21,*) energy
   energy = energy + mass * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT**2 / EL_CHG

   close(21)
else
   write(*,*) "ERROR, could not locate pncr file"
   stop
end if

ndata      = 2*norb*10     !< Maximum number of data points to be stored
fnout      = "poincare.h5" !< Output file

call random_seed()
call sim%initialize(num_groups=1)

if( mass .eq. 0) then
   allocate(particle_fieldline::sim%groups(1)%particles(nprt))
else
   allocate(particle_gc_relativistic::sim%groups(1)%particles(nprt))
end if

! Allocate and initialize needed data arrays
allocate ( rvals(nprt,ndata), zvals(nprt,ndata), phivals(nprt,ndata), psivals(nprt,ndata), &
     milvals(nprt, ndata), pncrid(nprt,ndata) )
allocate ( nextslot(nprt), ncross(nprt), mileage(nprt) )

call cpu_time(t0)

! Begin simulation at new time slice by initializing the field at that point
sim%time = 0
fieldreader = event(read_jorek_fields_interp_linear(i=-1))
call with(sim,fieldreader)

call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

! Find the location of the magnetic axis and initialize new markers
call find_axis(sim%my_id, sim%fields%node_list, sim%fields%element_list, &
     psiaxis, raxis, zaxis, ielmaxis, saxis, taxis, ifail)
call init_markers(sim, nprt, raxis, zaxis, mass, charge, pitch, energy)
write(*,*) "Markers initialized. Begin simulation."

! Init data arrays
rvals    = 0
zvals    = 0
psivals  = 0
phivals  = 0
pncrid   = 0
nextslot = 0
ncross   = 0
mileage  = 0

finished = 0
do while(finished .lt. nprt)

   finished = 0

   select type (p => sim%groups(1)%particles)
   type is (particle_gc_relativistic)

#ifdef __GFORTRAN__
      !$omp parallel do default(shared) & !To avoid GNU compiler failure
#else
      !$omp parallel do default(none) &
#endif
#ifdef __GFORTRAN__
      !$omp shared(sim, tstep, norb, ndata, mileage) & !To avoid GNU compiler failure
#else
      !$omp shared(sim, tstep, norb, ndata, mileage, p) &
#endif
      !$omp shared(raxis, zaxis, ncross, nprt, rvals, zvals, phivals, psivals, milvals, nextslot, pncrid) &
      !$omp private(iprt, ifail, zprev, phiprev, xyz, xyz0) &
      !$omp reduction(+:finished)
      do iprt=1,nprt
         if (p(iprt)%i_elm .gt. 0 .and. ncross(iprt) .lt. norb) then
            xyz0    = cylindrical_to_cartesian(p(iprt)%x)
            zprev   = p(iprt)%x(2)
            phiprev = p(iprt)%x(3)
            call runge_kutta_fixed_dt_gc_push_jorek(sim%fields,sim%time, tstep, &
                 sim%groups(1)%mass, p(iprt))

            xyz = cylindrical_to_cartesian(p(iprt)%x)
            if ( p(iprt)%i_elm > 0  ) then
              mileage(iprt) = mileage(iprt) + norm2(xyz - xyz0)
            endif
            call check_and_store_crossing(sim%fields, iprt, sim%time, mileage(iprt), p(iprt)%i_elm, p(iprt)%st, p(iprt)%x, &
                 phiprev, zprev, raxis, zaxis, ndata, nextslot, ncross, rvals, zvals, phivals, psivals, milvals, pncrid)
         else
            finished = finished + 1
         end if

      end do
      !$omp end parallel do
      
   type is (particle_fieldline)
      
#ifdef __GFORTRAN__
      !$omp parallel do default(shared) & ! To avoid GNU compiler failure
#else
      !$omp parallel do default(none) &
#endif
#ifdef __GFORTRAN__
      !$omp shared(sim, tstep, norb, ndata, mileage) &
#else
      !$omp shared(sim, tstep, norb, ndata, mileage, p) & ! To avoid GNU compiler failure
#endif
      !$omp shared(raxis, zaxis, ncross, nprt, rvals, zvals, phivals, psivals, milvals, nextslot, pncrid) &
      !$omp private(iprt, ifail, zprev, phiprev) &
      !$omp reduction(+:finished)
      do iprt=1,nprt
         if (p(iprt)%i_elm .gt. 0 .and. ncross(iprt) .lt. norb) then
            zprev   = p(iprt)%x(2)
            phiprev = p(iprt)%x(3)
            call field_line_runge_kutta_fixed_dt_push_jorek(sim%fields, p(iprt), sim%time, tstep)

            mileage(iprt) = mileage(iprt) + tstep
            call check_and_store_crossing(sim%fields, iprt, sim%time, mileage(iprt), p(iprt)%i_elm, p(iprt)%st, p(iprt)%x, &
                 phiprev, zprev, raxis, zaxis, ndata, nextslot, ncross, rvals, zvals, phivals, psivals, milvals, pncrid)
         else
            finished = finished + 1
         end if

      end do
      !$omp end parallel do

   end select

end do

call cpu_time(t1)
write(*,*) 'CPU time: ', t1-t0

! Store data
call write_poincare_hdf5(fnout, rvals, zvals, phivals, psivals, milvals, pncrid, mileage)

! Finalize the simulation
deallocate ( rvals, zvals, phivals, psivals, milvals, pncrid, nextslot, ncross, mileage )
call sim%finalize


contains

!< Initialize markers at the outer mid plane
!<
!< Markers are initialized in z = zaxis plane and distributed between [raxis, redge] in fixed
!< intervals. Toroidally markers are randomly uniformly distributed.
!<
!< The type of sim%groups%particles decides what markers are initialized (field lines or
!< guiding centers). For physical particles mass, pitch, and energy are required.
subroutine init_markers(sim, nprt, raxis, zaxis, mass, charge, pitch, energy)
  implicit none

  type(particle_sim), intent(inout) :: sim    !< The particle simulation struct
  integer*4, intent(in)             :: nprt   !< Numebr of markers to be initialized
  real*8, intent(in)                :: raxis  !< Axis R coordinate in m
  real*8, intent(in)                :: zaxis  !< Axis z coordinate in m
  real*8, intent(in)                :: mass   !< Particle mass in AMU
  integer*4, intent(in)             :: charge !< Particle charge in e
  real*8, intent(in)                :: pitch  !< Particle pitch (vpar/v)
  real*8, intent(in)                :: energy !< Particle energy in eV including the rest energy
  
  integer*4 :: iprt, ifail, i_elm
  real*8, allocatable :: rnd(:)   !< Array for storing random numbers
  real*8 :: redge, st(2)

  type(particle_gc_relativistic) :: particle_in, particle_out

  !! Find redge (within accuracy of 1 cm) !!
  redge = raxis
  ifail = 0
  do while(ifail .eq. 0)
     redge = redge + 0.1d0
     call find_RZ(sim%fields%node_list, sim%fields%element_list, &
          redge, zaxis, redge, zaxis, i_elm, st(1), st(2), ifail)
  end do
  redge = redge - 0.1d0


  !! Initialize markers !!
  sim%groups(1)%mass = mass

  allocate ( rnd(nprt) )
  call random_number(rnd)

  do iprt = 1,nprt
     select type (p=>sim%groups(1)%particles(iprt))
     type is (particle_fieldline)

        p%x = [raxis + iprt * ( redge - raxis ) / (nprt+1), zaxis, 2*PI*rnd(iprt)]

        call find_RZ(sim%fields%node_list, sim%fields%element_list, &
             p%x(1), p%x(2), &
             p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail)

        p%v = 1e0 ! field line "velocity" is set to 1 m/s.

     type is (particle_gc_relativistic)
        p%q = charge                                                                                                                                                                                                                                                                     
        p%x = [raxis + iprt * ( redge - raxis ) / (nprt+1), zaxis, 2*PI*rnd(iprt)]

        call find_RZ(sim%fields%node_list, sim%fields%element_list, &
             p%x(1), p%x(2), &
             p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail)

        if (p%i_elm .gt. 0) then
           particle_out = relativistic_gc_momenta_from_E_cospitch(p, energy, pitch, mass, sim%fields, sim%time)
           p%p = particle_out%p
        end if
     end select
  end do

  deallocate ( rnd )
   
end subroutine init_markers


!< Check whether Poincare plane was crossed during this time step and store crossing if this was the case.
!<
!< Note that this method actually records the position marker has *after* it has crossed the plane, i.e., the
!< actual location of the crossing is not stored. This method could be improved by performing a linear interpolation
!< on the marker positions at each side of the plane.
subroutine check_and_store_crossing(fields, iprt, time, mileage, i_elm, st, x, phiprev, zprev, raxis, zaxis, ndata, &
     nextslot, ncross, rvals, zvals, phivals, psivals, milvals, pncrid)
  implicit none
  class(fields_base), intent(in) :: fields !< The field structure for evaluating psi
  integer*4, intent(in) :: iprt    !< Marker position in the sim%group array
  real*8, intent(in)    :: time    !< Current time
  real*8, intent(in)    :: mileage !< Distance marker has travelled
  integer*4, intent(in) :: i_elm   !< Element the marker is located in
  real*8, intent(in)    :: st(2)   !< Particle position within the element
  real*8, intent(in)    :: x(3)    !< Marker R,z,phi coordinates
  real*8, intent(in)    :: phiprev !< Marker toroidal coordinate at previous time step
  real*8, intent(in)    :: zprev   !< Marker z coordinate at previous time step
  real*8, intent(in)    :: raxis   !< Magnetic axis R coordinate
  real*8, intent(in)    :: zaxis   !< Magnetic axis z coordinate
  integer*4, intent(in) :: ndata   !< Maximum number of data points
  
  integer*4, intent(inout) :: nextslot(:)  !< Array (len=nprt) storing the index of next available slot
  integer*4, intent(inout) :: ncross(:)    !< Array (len=nprt) storing the number of crossings
  real*8, intent(inout)    :: rvals(:,:)   !< Stored R coordinate of a crossing (nprt,ndata)
  real*8, intent(inout)    :: zvals(:,:)   !< Stored z coordinate of a crossing (nprt,ndata)
  real*8, intent(inout)    :: phivals(:,:) !< Stored phi coordinate of a crossing (nprt,ndata)
  real*8, intent(inout)    :: psivals(:,:) !< Stored psi_n coordinate of a crossing (nprt,ndata)
  real*8, intent(inout)    :: milvals(:,:) !< Stored mileage of a crossing (nprt,ndata)
  integer*4, intent(inout) :: pncrid(:,:)  !< Stored plane index of a crossing (nprt,ndata)

  real*8 :: E(3), B(3), U, psi, psin
  integer*4 :: idx !< Helper variable

  ! Check whether OMP was crossed and store crossing
  if (i_elm .gt. 0 .and. &
       ( x(1) > raxis .and. sign(1.0, x(2) - zaxis) .ne. sign(1.0, zprev - zaxis) ) ) then
     ncross(iprt)       = ncross(iprt) + 1 ! Update the turns counter

     nextslot(iprt)     = nextslot(iprt) + 1
     if(nextslot(iprt) .gt. ndata) then
        nextslot(iprt)  = 1
     end if

     call get_psi_n0(fields, i_elm, st(1), st(2), x(2), psin)
     
     idx = nextslot(iprt)
     rvals(iprt, idx)   = x(1)
     zvals(iprt, idx)   = x(2)
     phivals(iprt, idx) = x(3)
     psivals(iprt, idx) = psin
     milvals(iprt, idx) = mileage
     pncrid(iprt, idx)  = 1
  end if

  ! Check whether the poloidal plane at phi = 0 was crossed and store the crossing
  if (i_elm .gt. 0 .and. floor(phiprev/(2*PI)) .ne. floor(x(3)/(2*PI)) ) then
     nextslot(iprt)     = nextslot(iprt) + 1
     if(nextslot(iprt) .gt. ndata) then
        nextslot(iprt)  = 1
     end if

     call get_psi_n0(fields, i_elm, st(1), st(2), x(2), psin)

     idx = nextslot(iprt)
     rvals(iprt, idx)   = x(1)
     zvals(iprt, idx)   = x(2)
     phivals(iprt, idx) = x(3)
     psivals(iprt, idx) = psin
     milvals(iprt, idx) = mileage
     pncrid(iprt, idx)  = 2
  end if

end subroutine check_and_store_crossing


!< Write the Poincare output to a HDF5 file
!<
!< The output is stored in 1D arrays. To dissect the data, one of the arrays (iprt)
!< identifies to which marker this data point corresponds to while one array (pncrid)
!< identifies the plane that was crossed (1 : outer mid-plane, 2: poloidal plane).
!<
!< The recorded points contain information on the cylindirical coordinates on each crossing.
!< There is an additional mileage array where mileage(iprt) notes the total distance iprt
!< marker travelled during the simulation (seconds for particles and meters for field lines).
!< Therefore, for lost markers mileage is the connection length.
subroutine write_poincare_hdf5(fnout, rvals, zvals, phivals, psivals, milvals, pncrid, mileage)
  implicit none
  character(len=40)     :: fnout        !< File name where the data is written e.g. "poincare.h5"
  real*8, intent(in)    :: rvals(:,:)   !< Recorded R coordinates (nprt, ndata)
  real*8, intent(in)    :: zvals(:,:)   !< Recorded z coordinates (nprt, ndata)
  real*8, intent(in)    :: phivals(:,:) !< Recorded phi coordinates (nprt, ndata)
  real*8, intent(in)    :: psivals(:,:) !< Recorded psi_n coordinates (nprt, ndata)
  real*8, intent(in)    :: milvals(:,:) !< Recorded mileage at each crossing (nprt, ndata)
  integer*4, intent(in) :: pncrid(:,:)  !< Plane ID for each crosssing (nprt, ndata)
  real*8, intent(in)    :: mileage(:)   !< Travelled distance (nprt)

  
  real*8, allocatable :: rvals0(:), zvals0(:)     !< 1D arrays for storing the data once empty points are removed
  real*8, allocatable :: phivals0(:), psivals0(:) !< -,,-
  real*8, allocatable :: milvals0(:)              !< -,,-
  integer*4, allocatable :: pncrid0(:), iprt(:)   !< Arrays for identifyinf the marker and plane from the 1D data arrays
  integer*4 :: nprt, ndata, npoint, i, j, ipoint  !< Helper variables
  integer(HID_T) :: file
  integer :: hdferr

  nprt  = size(rvals,1)
  ndata = size(rvals,2)

  ! Find out how many points with actual data there is
  npoint = 0
  do i=1,nprt
     do j=1,ndata
        if(pncrid(i,j) .ne. 0) then
           npoint = npoint + 1
        end if
     end do
  end do

  allocate( iprt(npoint), pncrid0(npoint) )
  allocate( rvals0(npoint), zvals0(npoint), phivals0(npoint), psivals0(npoint), milvals0(npoint) )

  ! Store the (actual) data to 1D arrays
  ipoint = 1
  do i=1,nprt
     do j=1,ndata
        if(pncrid(i,j) .ne. 0) then
           iprt(ipoint)     = i
           pncrid0(ipoint)  = pncrid(i,j)
           rvals0(ipoint)   = rvals(i,j)
           zvals0(ipoint)   = zvals(i,j)
           phivals0(ipoint) = mod(phivals(i,j), 2*PI)
           psivals0(ipoint) = psivals(i,j)
           milvals0(ipoint) = milvals(i,j)
           ipoint = ipoint + 1
        end if
     end do
  end do

  ! Write to HDF5
  call h5open_f(hdferr)
  call h5fcreate_f(fnout, H5F_ACC_TRUNC_F, file, hdferr)
  if (hdferr .gt. 0) then
     write(*,*) "file open failed:", hdferr
     return
  end if

  call HDF5_array1D_saving(file, rvals0,   npoint, "r")
  call HDF5_array1D_saving(file, zvals0,   npoint, "z")
  call HDF5_array1D_saving(file, phivals0, npoint, "phi")
  call HDF5_array1D_saving(file, psivals0, npoint, "psi")
  call HDF5_array1D_saving(file, milvals0, npoint, "mil")
  call HDF5_array1D_saving_int(file, pncrid0,  npoint, "pncrid")
  call HDF5_array1D_saving_int(file,    iprt,  npoint, "iprt")
  call HDF5_array1D_saving(file, mileage, nprt, "mileage")

  deallocate( iprt, pncrid0, rvals0, zvals0, phivals0, psivals0, milvals0 )

  ! This would be the raw data that is in 2D arrays that someone might find easier to handle
  ! and therefore it is left here
  !call HDF5_array2D_saving(file, rvals,   nprt, norb, "r")
  !call HDF5_array2D_saving(file, zvals,   nprt, norb, "z")
  !call HDF5_array2D_saving(file, phivals, nprt, norb, "phi")
  !call HDF5_array2D_saving(file, psivals, nprt, norb, "psi")
  !call HDF5_array2D_saving(file, milvals, nprt, norb, "mil")
  !call HDF5_array2D_saving_int(file, pncrid,  nprt, norb, "pncrid")

  call h5fclose_f(file, hdferr)
  call h5close_f(hdferr)

end subroutine write_poincare_hdf5


!< Helper function to evaluate normalized psi from n = 1 component
subroutine get_psi_n0(field, i_elm, s, t, z, psin)
  implicit none

  class(fields_base),  intent(in) :: field !< Field data
  integer, intent(in) :: i_elm             !< Corresponding element
  real*8, intent(in)  :: s, t, z           !< Position in element coordinates and z
  real*8, intent(out) :: psin              !< Normalized psi

  real*8 :: P(1), P_s, P_t, P_st, P_ss, P_tt
  integer :: i_var = 1, n_tor = 1

  call interp(field%node_list, field%element_list, i_elm, i_var, n_tor, s, t, P(1), P_s, P_t, P_st, P_ss, P_tt)
  psin = get_psi_n(P(1), z)

end subroutine get_psi_n0

end program poincare
