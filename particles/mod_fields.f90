!> Module containing base type for field interpolations, with interfaces
!> to implement
module mod_fields
use data_structure
implicit none
private
public fields_base

!> Base type for a field interpolator.
!> Must implement the following interfaces, which are the normal
!> functions and an additional time component (JOREK units)
!> node_list and element_list should be the currently-valid representation of the grid
!> (values themselves should not be used, only for find_RZ etc)
type, abstract :: fields_base
  type(type_node_list),pointer         :: node_list => null()!< Current node list
  type(type_element_list), pointer     :: element_list => null() !< Current element list
  logical                              :: static=.false. !< if true do not time interpolate
  logical                              :: flag_zero_dpsidt=.false. !< if true, P_time(1) = dpsi/dt = 0
  contains
    procedure(interp_PRZ), deferred, public   :: interp_PRZ
    procedure(interp_PRZ_2), deferred, public :: interp_PRZ_2
    procedure, public :: calc_NeTe
    procedure, public :: calc_EBpsiU
    procedure, public :: calc_EBNormBGradBCurlbDbdt
    procedure, public :: calc_analytical_EBpsiU
    procedure, public :: calc_analytical_EBNormBGradBCurlbDbdt
    procedure, public :: set_flag_dpsidt
end type fields_base

interface
  !> Interpolate a variable at s, t, phi in i_elm, returning first
  !> derivatives of the variable and of space
  pure subroutine interp_PRZ(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
    import fields_base
    class(fields_base),  intent(in)  :: this
    real*8,                   intent(in)  :: time !< Time at which to calculate this variable
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_time(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    real*8,                   intent(out) :: P_phi(n_v)
  end subroutine interp_PRZ
  !> Interpolate a variable at s, t, phi in i_elm, returning first
  !> and second order derivatives of the variable and of R and Z.
  pure subroutine interp_PRZ_2(this,time,i_elm,i_v,n_v,s,t,phi,P,P_s,P_t,P_phi, &
    P_time,P_ss,P_st,P_tt,P_sphi,P_tphi,P_stime,P_ttime,                        &
    R,R_s,R_t,R_ss,R_st,R_tt,Z,Z_s,Z_t,Z_ss,Z_st,Z_tt)
    import fields_base
    !> declare input variables
    class(fields_base), intent(in)      :: this
    real(kind=8), intent(in)            :: time, s, t, phi
    integer, intent(in)                 :: i_elm, n_v
    integer, dimension(n_v), intent(in) :: i_v
    !> declare ourput variables
    real(kind=8), intent(out)                 :: R, R_s, R_t, R_ss, R_st, R_tt
    real(kind=8), intent(out)                 :: Z, Z_s, Z_t, Z_ss, Z_st, Z_tt
    real(kind=8), dimension(n_v), intent(out) :: P, P_s, P_t, P_phi, P_time
    real(kind=8), dimension(n_v), intent(out) :: P_ss, P_st, P_tt, P_sphi, P_tphi
    real(kind=8), dimension(n_v), intent(out) :: P_stime, P_ttime
  end subroutine interp_PRZ_2
end interface

contains
!> Calculates the electric and magnetic fields at a specific position
!> in the jorek element `i_elm` at `st`.
pure subroutine calc_EBpsiU(fields, time, i_elm, st, phi, E, B, psi, U)
use phys_module, only: F0, mode, central_mass, central_density
use constants, only: mu_zero, mass_proton
use mod_coordinate_transforms, only: transform_derivatives_st_to_RZ
! Routine parameters
class(fields_base), intent(in) :: fields
real*8, intent(in)  :: time
integer, intent(in) :: i_elm !< JOREK element index
real*8, intent(in)  :: st(2) !< element-local coordinates
real*8, intent(in)  :: phi !< toroidal angle
real*8, intent(out) :: E(3) !< Electric field [V/m]
real*8, intent(out) :: B(3) !< Magnetic field [T]
real*8, intent(out) :: psi !< psi in JOREK units
real*8, intent(out) :: u !< velocity stream function in m/s

! Internal parameters
integer, parameter :: i_var(2) = [1,2]
real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_time(2) ! Placeholder for evaluating variables and derivatives locally
! Values
real*8             :: R, R_s, R_t, Z, Z_s, Z_t
! Others
real*8             :: inv_st_jac, R_inv
real*8             :: psi_R, psi_Z, U_R, U_Z, U_phi, t_norm

t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

! Interpolate the fields to get psi and U at the current position (and the
! changes u_n - u(n-1))
call fields%interp_PRZ(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

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
B     = [+psi_Z, -psi_R, F0] * R_inv

! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
! See http://jorek.eu/wiki/doku.php?id=u_phi
E     = [-F0*U_R, -F0*U_Z, -F0*U_phi*R_inv]/t_norm
E(3)  = E(3) - R_inv*P_time(1) ! because this is not normalized with t_norm

end subroutine calc_EBpsiU

pure subroutine calc_NeTe(fields, time, i_elm, st, phi, n_e, T_e, grad_T_e)
use phys_module, only: central_density
use constants
class(fields_base), intent(in)                    :: fields
integer, intent(in)                               :: i_elm
real*8, intent(in)                                :: time, st(2), phi
real*8, intent(out)                               :: n_e !< electron density [m^-3]
real*8, intent(out)                               :: T_e !< electron temperature [K]
real*8, intent(out), optional, dimension(3)       :: grad_T_e !< gradient of electron temperature [K/m]

real*8, dimension(2) :: P, P_s, P_t, P_phi, P_time
real*8               :: R, R_s, R_t, Z, Z_s, Z_t, xjac
real*8 :: T_norm !< temperature normalisation

call fields%interp_PRZ(time,i_elm,&
#if (JOREK_MODEL == 400)
      [5,8],& ! electron temperature
#else
      [5,6],& ! electron temperature + ion temperature (assumed equal)
#endif
          2,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)

n_e = max(central_density * P(1) * 1d20,1d16)                           ! plasma density [1/m^3], capped against negative
T_norm = (1.d0/K_BOLTZ/(2.d0*MU_ZERO*central_density*1.d20))
#if (JOREK_MODEL == 400)
T_norm = T_norm*2.d0 ! P(1) contains the electron temperature, reverse previous correction
#endif
T_e = max(P(2)*T_norm, 1.d0) ! temperature capped against going negative

if (present(grad_T_e)) then

  xjac = R_s * Z_t - R_t * Z_s
  grad_T_e = T_norm*[(  P_s(2) * Z_t - P_t(2) * Z_s)/ xjac, &
                     (- P_s(2) * R_t + P_t(2) * R_s)/ xjac, &
                     P_phi(2)/R]
end if
end subroutine calc_NeTe

!> This procedure computes the fields appearing in the
!> the guiding center equations of motion
!> inputs:
!>   fields: (field_base) structure containing methods for computing EM fields
!>   time:   (real8) current time
!>   i_elm:  (integer) particle mesh element index
!>   st:     (real8) particle position in local mesh coordinates
!>   phi:    (real8) particle toroidal angle
!> outputs:
!>   E:      (real8)(3) electric field in V/m
!>   b:      (real8)(3) magnetic field direction
!>   normB:  (real8)(3) magnetic field intensity in T
!>   gradB:  (real8)(3) gradient of the magnetic field intensity in T/m
!>   curlb:  (real8)(3) curl of the magnetic field direction in 1/m
!>   dbdt:   (real8)(3) magnetic field direction time derivative in 1/s
pure subroutine calc_EBNormBGradBCurlbDbdt(fields,time,i_elm,st,phi,E,b, &
  normB,gradB,curlb,dbdt)
  !> load modules
  use phys_module, only: F0, mode, central_mass, central_density
  use constants, only: mu_zero,mass_proton
  use mod_math_operators, only: cross_product 
  use mod_coordinate_transforms, only: transform_first_derivatives_st_to_RZ
  use mod_coordinate_transforms, only: transform_second_derivatives_st_to_RZ
  implicit none

  !> declare input variables
  class(fields_base), intent(in)         :: fields
  real(kind=8), intent(in)               :: time
  integer, intent(in)                    :: i_elm
  real(kind=8), dimension(2), intent(in) :: st
  real(kind=8), intent(in)               :: phi
  !> declare output variables
  real(kind=8), intent(out)               :: normB
  real(kind=8), dimension(3), intent(out) :: E, b, gradB, curlb, dbdt
  !> declare internal variables
  real(kind=8)               :: R_inv, normB_inv
  real(kind=8), dimension(2) :: U_RZ !< 1:U_R, 2:U_Z
  !> global coordinates and derivatives: 1:R, 2:R_s, 3:R_t, 4:R_ss, 5:R_st, 6:R_tt,
  !> 7:Z, 8:Z_s, 9:Z_t, 10:Z_ss, 11:Z_st, 12:Z_tt
  real(kind=8), dimension(12) :: RZ
  real(kind=8), dimension(5)  :: U !< stream function: [U,U_R,U_Z,U_phi,U_time]
  !> psi derivatives in global coordinates: 1:psi_R, 2:psi_Z, 3:psi_RR,
  !> 4:psi_RZ, 5:psi_ZZ, 6:psi_Rphi, 7:psi_Zphi, 8:psi_Rtime, 9:psi_Ztime
  real(kind=8), dimension(9) :: psi_RZ
  !> poloidal flux: psi, psi_s, psi_t, psi_phi, psi_time, psi_ss, psi_st, psi_tt,
  !>   psi_sphi, psi_tphi, psi_stime, psi_ttime
  real(kind=8), dimension(12) :: psi 

  !> interpolate the stream function
  call fields%interp_PRZ(time,i_elm,[2],1,st(1),st(2),phi,U(1),U(2),U(3), &
    U(4),U(5),RZ(1),RZ(2),RZ(3),RZ(7),RZ(8),RZ(9))

  !> convert the electric potential into SI units
  U = F0*U/sqrt(mu_zero*mass_proton*central_mass*central_density*1.d20)

  !> transform first U derivatives from st to RZ
  call transform_first_derivatives_st_to_RZ(U_RZ(1),U_RZ(2),1,U(2),U(3), &
    RZ(2),RZ(3),RZ(8),RZ(9))
  
  !> interpolate the poloidal flux
  call fields%interp_PRZ_2(time,i_elm,[1],1,st(1),st(2),phi,psi(1),psi(2),&
       psi(3),psi(4),psi(5),psi(6),psi(7),psi(8),psi(9),psi(10),psi(11),&
       psi(12),RZ(1),RZ(2),RZ(3),RZ(4),RZ(5),RZ(6),RZ(7),RZ(8),RZ(9),&
       RZ(10),RZ(11),RZ(12))

  !> set dpsidt to zero if needed
  if(fields%flag_zero_dpsidt) then
     psi(5)  = 0.d0 !< psi_time
     psi(11) = 0.d0 !< psi_stime
     psi(12) = 0.d0 !< psi_ttime
  endif

  R_inv = 1.d0/RZ(1) !< compute the inverse of R

  !> transform first and second order psi derivatives from st to RZ
  call transform_first_derivatives_st_to_RZ(psi_RZ(1),psi_RZ(2),1,psi(2),psi(3),&
       RZ(2),RZ(3),RZ(8),RZ(9))
  call transform_second_derivatives_st_to_RZ(psi_RZ(3),psi_RZ(4),psi_RZ(5),1,&
       psi(6),psi(7),psi(8),psi_RZ(1),psi_RZ(2),RZ(2),RZ(3),RZ(4),RZ(5),RZ(6),&
       RZ(8),RZ(9),RZ(10),RZ(11),RZ(12))
  call transform_first_derivatives_st_to_RZ(psi_RZ(6),psi_RZ(7),1,psi(9),psi(10),&
       RZ(2),RZ(3),RZ(8),RZ(9))
  call transform_first_derivatives_st_to_RZ(psi_RZ(8),psi_RZ(9),1,psi(11),psi(12),&
       RZ(2),RZ(3),RZ(8),RZ(9))
  
  !> compute the electric field
  E = -[U_RZ(1),U_RZ(2),R_inv*(U(4)+psi(5))] !< V/m

  !> compute the magnetic field (put it in variable b temporarily to save one variable)
  b = [psi_RZ(2),-psi_RZ(1),F0]*R_inv !< magnetic field T

  normB = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3)) !< B field intensity

  normB_inv = 1.d0/normB !< inverse of the B field intensity

  !< direction of the magnetic field
  b = b/normB

  !> compute the gradB field
  gradB = [psi_RZ(1)*psi_RZ(3)+psi_RZ(2)*psi_RZ(4),                          &
    psi_RZ(1)*psi_RZ(4)+psi_RZ(2)*psi_RZ(5),                                 &
    R_inv*(psi_RZ(1)*psi_RZ(6)+psi_RZ(2)*psi_RZ(7))]*R_inv*R_inv*normB_inv
  gradB(1) = gradB(1)-normB*R_inv
  
  !> compute the curlb field
  curlb = normB_inv*(cross_product(b,gradB) + &
    R_inv*[R_inv*psi_RZ(6),R_inv*psi_RZ(7),   &
    R_inv*psi_RZ(1)-psi_RZ(3)-psi_RZ(5)])

  !> compute the dbdt field
  dbdt = ((b(2)*psi_RZ(8)-b(1)*psi_RZ(9))*b +    &
    [psi_RZ(9),-psi_RZ(8),0.d0])*normB_inv*R_inv
  
end subroutine calc_EBNormBGradBCurlbDbdt

!> Subroutine to ocompute analytical magnetic and electric fields
!> for testing integrators. The electric field is set to zero
!> while a tokamak-like magnetic field with a poloidal flux of
!> 0.5*B0*((R-R0)**2+(Z-Z0)**2) is used.
!> inputs:
!>   RZ: (real8) particle poloidal plane position
!> outputs:
!>   B:   (real8)(3) magnetic field
!>   E:   (real8)(3) electric field
!>   psi: (real8) poloidal flux
pure subroutine calc_analytical_EBpsiU(fields,RZ,E,B,psi,U)
  implicit none
  !> declare parameters
  real(kind=8), parameter :: B0=2.5d0 !< axis magnetic field in [T]
  real(kind=8), parameter :: U0=0.d0 !< reference electric potential
  !> set magnetic axis position
  real(kind=8), dimension(2), parameter :: RZ0=[3.d0,0.d0]
  !> delcare input variables
  class(fields_base), intent(in) :: fields
  real(kind=8), dimension(2), intent(in) :: RZ
  !> declare output variables:
  real(kind=8), intent(out) :: psi, U
  real(kind=8), dimension(3), intent(out) :: E, B

  !> computing magnetic field
  B = B0*[RZ(2)-RZ0(2),RZ0(1)-RZ(1),RZ0(1)]/RZ(1)

  !> computing electric field
  E = U0*[0.d0,0.d0,0.d0]

  !> compute psi
  psi = 0.5*B0*(dot_product(RZ-RZ0,RZ-RZ0))

  !> compute U
  U = U0
  
end subroutine calc_analytical_EBpsiU

!> This procedure computes analytical guiding ceneter
!> fields for a static electromagnetic field. The
!> electric field is set to zero while a tokamak-like
!> magnetic field with a poloidal flux of:
!> psi = 0.5*B0*((R-R0)**2 + (Z-Z0)**2) is used.
!> inputs:
!>   RZ: (real8)(2) particle position in the poloidal plane
!> outputs:
!>   E:     (real8)(3) electric field
!>   b:     (real8)(3) magnetic field direction
!>   normB: (real8) magnetic intensity
!>   gradB: (real8)(3) gradient of the magnetic intensity
!>   curlb: (real8)(3) curl of the magnetic direction
!>   dbdt:  (real8)(3) magnetic direction time variation
pure subroutine calc_analytical_EBNormBGradBCurlbDbdt(fields, &
  RZ,E,b,normB,gradB,curlb,dbdt)
  use mod_math_operators, only: cross_product
  implicit none
  !> define parameters
  real(kind=8), parameter               :: B0=2.5d0 !< axis magnetic field in [T]
  real(kind=8), parameter               :: U0=0.d0  !< reference electric potential
  real(kind=8), dimension(2), parameter :: RZ0=[3.d0,0.d0]
  !> input variables
  class(fields_base), intent(in)         :: fields
  real(kind=8), dimension(2), intent(in) :: RZ
  !> output variables
  real(kind=8), intent(out)               :: normB
  real(kind=8), dimension(3), intent(out) :: E, b, gradB, curlb, dbdt

  !> compute electric field
  E = U0*[0.d0,0.d0,0.d0]
  
  !> compute magnetic field
  b = B0*[RZ(2)-RZ0(2),RZ0(1)-RZ(1),RZ0(1)]/RZ(1)
  
  !> compute norm of the magnetic field
  normB = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
  
  !> compute gradient of the magnetic field
  gradB = [B0*B0*(RZ(1)-RZ0(1))-normB*normB*RZ(1), &
    B0*B0*(RZ(2)-RZ0(2)),0.d0]/(normB*RZ(1)*RZ(1))

  !> compute the magetic direction
  b = b/normB

  !> compute the curl of the magnetic field directon
  curlb = (cross_product(b,gradB) -                 &
    [0.d0,0.d0,(RZ(1)+RZ0(1))/(RZ(1)*RZ(1))])/normB

  !> compute magnetic field time derivative
  dbdt = [0.d0,0.d0,0.d0]
  
end subroutine calc_analytical_EBNormBGradBCurlbDbdt

! This subroutine sets a flag to force dpsi/dt to 0
pure subroutine set_flag_dpsidt(this,flag_dpsidt_to_zero)
  class(fields_base),intent(inout) :: this !< fields object
  logical,intent(in)               :: flag_dpsidt_to_zero !< flag value

  this%flag_zero_dpsidt = flag_dpsidt_to_zero
  
end subroutine set_flag_dpsidt

end module mod_fields
