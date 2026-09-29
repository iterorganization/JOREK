!> Calculate the electric and magnetic field over time, from the two main
!> contributions of dA/dt and grad_// u.
!>
!> ./E_diagnostic start_time tstep psi_n < input
!>
!> Take a starting time and a timestep on the commandline, and timestep a
!> static particle set through it. These particles are located on a specific
!> flux surface given on the commandline. Uniformly distributed in phi
!> and theta.
!>
!> No MPI support
module write_E_time
  use mod_event, only: event, action
  use equil_info, only:find_xpoint
  use mod_interp
  implicit none

  type, extends(action) :: save_E
    character(len=20) :: suffix = '.dat'
  contains
    procedure :: do => do_save_E
  end type save_E

contains
!> Write out the contributions to n_var files for the particle positions in group 1
subroutine do_save_E(this, sim, ev)
  use mod_particle_sim, only: particle_sim
  use phys_module, only: central_mass, central_density, F0, eta, xpoint, xcase
  use constants, only: mu_zero, atomic_mass_unit
  class(save_E), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev

  integer, parameter :: n_var = 14
  character(len=7) :: variables(n_var) = ["E_r    ", "E_z    ", "E_phi  ", "E_phi2 ", &
                                          "B_r    ", "B_z    ", "B_phi  ", &
                                          "J_r    ", "J_z    ", "J_phi  ", "J_0    ", &
                                          "E_par  ", "E_par2 ", "E_par3 "]
  integer :: i, j, k, i_elm, i_elm_xpoint(2), ifail
  integer :: unit(n_var)
  real*8, dimension(4) :: P, P_s, P_t, P_phi, P_time
  real*8 :: R, R_s, R_t, Z, Z_s, Z_t
  real*8 :: psi_R, psi_Z, U_R, U_Z, U_phi
  real*8 :: R_inv, inv_st_jac
  real*8 :: t_norm
  real*8 :: E(3), B(3), E_phi
  real*8 :: J_0, J_phi, J_r, J_z, eta_T, T_ref
  real*8 :: s, t, s_xpoint(2), t_xpoint(2)
  real*8 :: psi_axis, psi_xpoint(2)
  real*8 :: R_axis, Z_axis, R_xpoint(2), Z_xpoint(2)
  character(len=12) time_s

  write(time_s,'(f12.8)') sim%time
  do i=1,n_var
    open(newunit=unit(i), file=trim(adjustl(time_s))//'_'//trim(variables(i))//this%suffix, status='replace', access='stream', action='write')
  end do

  t_norm  = sqrt(mu_zero * ATOMIC_MASS_UNIT * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  ! calculate current J0 to compare jphi - jphi0, since that is what we solve for
  call find_axis(0, sim%fields%node_list, sim%fields%element_list, psi_axis, R_axis, Z_axis, i_elm, s, t, ifail)
  ! T_ref = T_axis
  call sim%fields%interp_PRZ(sim%time, i_elm, [6], 1, &
    s, t, 0.d0, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
  T_ref = P(1)

  call find_xpoint(0, sim%fields%node_list, sim%fields%element_list, psi_xpoint, R_xpoint, Z_xpoint, i_elm_xpoint, s_xpoint, t_xpoint, xcase, ifail)

  do i=1,1 ! do not support multiple groups for now
    do j=1,size(sim%groups(i)%particles,1)
      if (sim%groups(i)%particles(j)%i_elm .gt. 0) then
        ! Interpolate the fields to get psi and U at the current position (and the
        ! changes u_n - u(n-1))
        call sim%fields%interp_PRZ(sim%time, sim%groups(i)%particles(j)%i_elm, [1,2,3,6], 4, &
          sim%groups(i)%particles(j)%st(1), sim%groups(i)%particles(j)%st(2), &
          sim%groups(i)%particles(j)%x(3), &
          P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

        R_inv = 1.d0/R
        inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)

        ! Calculate the derivatives to R and Z
        psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
        psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
        U_R      = (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
        U_Z      = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
        U_phi    = P_phi(2)

        ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
        ! B = F0 a^3 + 1/R dpsi/dz e_r - 1/R dpsi_dr e_z
        B     = [+psi_Z, -psi_R, F0] * R_inv

        ! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
        ! See http://jorek.eu/wiki/doku.php?id=u_phi
        E     = [-F0*U_R, -F0*U_Z, -F0*U_phi*R_inv]/t_norm
        E(3)  = E(3) - R_inv*P_time(1) ! separate because this is not normalized with t_norm

        ! use psi_xpoint(1) as psi_bnd here, so does not work for other xcases?
        call current(xpoint,xcase,R,Z,Z_xpoint(1),P(1),psi_axis,psi_xpoint(1),J_0)
        J_0 = -J_0*R_inv/MU_ZERO
        ! Total toroidal current Jphi - J0 
        ! J_phi in SI units = -J_jorek / R
        J_phi = (-P(3))*R_inv/MU_ZERO - J_0

        ! Poloidal current derived from psi
        ! F is assumed constant, so disregard the derivatives here
        ! J_r = 1/R^2 d_R d_phi psi
        ! J_z = 1/R^2 d_Z d_phi psi
        ! See https://www.jorek.eu/wiki/doku.php?id=reduced_mhd#current
        ! This is not available, since mod_fields_hermite_birkhoff does not support second derivatives yet
        J_R = 0.d0!R_inv*R_inv/MU_ZERO * 
        J_Z = 0.d0!R_inv*R_inv/MU_ZERO * 


        eta_T = eta * (T_ref/P(4))*sqrt(T_ref/P(4)) * MU_ZERO/t_norm

        ! v x B - eta J
        ! todo check signs
        E_phi = -(psi_R * U_Z - psi_Z * U_R)/t_norm &
                + eta_T * J_phi
        ! Compare E_phi here with E_phi from eta_J and poisson bracket

        write(unit(1)) real(E(1),4)
        write(unit(2)) real(E(2),4)
        write(unit(3)) real(E(3),4)
        write(unit(4)) real(E_phi,4) ! = E_phi2, the alternative calculation

        write(unit(5)) real(B(1),4)
        write(unit(6)) real(B(2),4)
        write(unit(7)) real(B(3),4)

        write(unit(8)) real(J_r,4)
        write(unit(9)) real(J_z,4)
        write(unit(10)) real(J_phi,4)
        write(unit(11)) real(J_0,4)

        write(unit(12)) real(dot_product(E,B)/norm2(B),4)
        write(unit(13)) real(dot_product([E(1), E(2), E_phi],B)/norm2(B),4)
        write(unit(14)) real(dot_product(eta_T*[J_r,   J_z, J_phi],B)/norm2(B),4)
        ! units for all are V/m, T or A/m^2
      else
        do k=1,n_var
          write(unit(k)) 0.0
        end do
      end if
    end do
  end do

  write(*,*) 'output written to ', trim(time_s), '_$variable.dat ', variables
  do i=1,n_var
    close(unit(i))
  end do

end subroutine do_save_E
end module write_E_time

program E_par_time
use write_E_time
use mod_particle_sim
use mod_fields_linear, only: last_file_before_time
use mod_fields_hermite_birkhoff
use phys_module, only: xpoint, xcase
use mod_random_seed, only: random_seed
use constants, only: TWOPI
use mod_event, only: action, event, next_event_at, with
use domains
use mod_particle_types, only: particle_fieldline
!$ use omp_lib
implicit none

integer, parameter :: n_theta = 2000, n_phi=1
! theta from 0 to 2pi (outboard midplane = 0, counter-clockwise)
! phi from 0 to 2pi

type(particle_sim) :: sim
type(event) :: fieldreader
integer :: i_elm, i_elm_xpoint(2), ifail, i, j, k, ierr
real*8 :: psi_n, delta_t
real*8 :: psi_axis, psi_xpoint(2)
real*8 :: R_axis, Z_axis, R_xpoint(2), Z_xpoint(2)
real*8, allocatable :: psi_minmax_list(:,:)
real*8 :: s, s_xpoint(2), t, t_xpoint(2), R, Z, phi, theta
character(len=20) :: time_s, suffix
type(event), allocatable :: events(:)

! Start up MPI, jorek
call sim%initialize(num_groups=1, skip_jorek2help=.true.)


call get_command_argument(1, time_s)
read(time_s,*) sim%time

call get_command_argument(2, time_s)
read(time_s,*) delta_t

call get_command_argument(3, time_s)
read(time_s,*) psi_n

if (psi_n .gt. 1 .or. psi_n .lt. 0) then
  write(*,*) 'psi_n input out of bounds', psi_n
  stop 1
end if

write(*,*) "Sampling from 0 to 2pi in phi (n=", n_phi, ") and then 0 to 2pi in theta (n=", n_theta, ") at psi_n = ", psi_n
write(*,*) "From t_start=", sim%time, " every dt=", delta_t, " s"


write(suffix,'(A,i5.5,A,i5.5,A)') '_', n_theta, '_', n_phi, '.dat'
! Use our event loop to run the diagnostic every delta_t
events = [event(read_jorek_fields_interp_hermite_birkhoff(i=last_file_before_time(sim%time))), &
          event(save_E(suffix=suffix), step=delta_t)]
call with(sim, events(1))



! Calculate normalization coefficients psi_axis and psi_xpoint(1)
! hack for asdex upgrade case since find_axis fails there
call find_axis(0, sim%fields%node_list, sim%fields%element_list, psi_axis, R_axis, Z_axis, i_elm, s, t, ifail)
call find_xpoint(0, sim%fields%node_list, sim%fields%element_list, psi_xpoint, R_xpoint, Z_xpoint, i_elm_xpoint, s_xpoint, t_xpoint, xcase, ifail)
write(*,*) "Psi_xpoint: ", psi_xpoint(1), " Psi_axis: ", psi_axis
if (abs(psi_xpoint(1) - psi_axis) .le. 1d-3) then
  write(*,*) "error in find_axis! same as xpoint... exiting"
  stop 1
end if

if (psi_xpoint(1) .eq. 0.d0 .or. psi_axis .eq. 0.d0) then
  write(*,*) "error in calculation of psi_xpoint or psi_axis"
  stop 1
end if




! Precalculate the positions
allocate(particle_fieldline::sim%groups(1)%particles(n_theta*n_phi))

allocate(psi_minmax_list(sim%fields%element_list%n_elements,2))
! Preparatory work: determine psi_min,max
!$omp parallel do default(shared) &
!$omp private(i_elm)
do i_elm=1,sim%fields%element_list%n_elements
  call psi_minmax(sim%fields%node_list,sim%fields%element_list,i_elm,psi_minmax_list(i_elm,1),psi_minmax_list(i_elm,2))
end do
!$omp end parallel do

open(unit=23, file='position.dat', access='stream', status='replace')
do i=1,n_theta
  do j=1,n_phi
    phi   = real(j-1,8)/real(n_phi,8) * TWOPI
    theta = real(i-1,8)/real(n_theta,8) * TWOPI
    call find_theta_psi(sim%fields%node_list,sim%fields%element_list,psi_minmax_list,&
      theta, psi_n*(psi_xpoint(1)-psi_axis) + psi_axis,phi,R_axis,Z_axis,i_elm,s,t,R,Z, ierr)
    if (ierr .ne. 0) then
      write(*,*) 'Problem finding position, setting 0'
      i_elm = 0
      s = 0.d0
      t = 0.d0
      R = 0.d0
      Z = 0.d0
    end if
    k = j + (i-1)*n_phi
    sim%groups(1)%particles(k)%i_elm = i_elm
    sim%groups(1)%particles(k)%st = [s,t]
    sim%groups(1)%particles(k)%x = [R,Z,phi]
    write(23) real(R,4), real(Z,4), real(theta,4), real(phi,4), real(i_elm,4), real(s,4), real(t,4)
  end do
end do
close(23)


do while (.not. sim%stop_now)
  sim%time = next_event_at(sim, events)
  call with(sim, events, at=sim%time)
end do
call sim%finalize
end program E_par_time
