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
  type(type_node_list),pointer         :: node_list    => null() !< Current node list
  type(type_element_list), pointer     :: element_list => null() !< Current element list
  logical                              :: static=.false. !< if true do not time interpolate
  logical                              :: flag_zero_dpsidt=.false. !< if true, P_time(1) = dpsi/dt = 0
  contains
    procedure(interp_PRZ), deferred, public   :: interp_PRZ
    procedure(interp_PRZ_2), deferred, public :: interp_PRZ_2
    procedure, public :: calc_NeTe
    procedure, public :: calc_NeTe_imp
    procedure, public :: calc_EBpsiU
    procedure, public :: calc_VBpsiU
    procedure, public :: calc_Qin, calc_Qin_analytic
    procedure, public :: calc_rk4, calc_RK4_analytic
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
subroutine calc_EBpsiU(fields, time, i_elm, st, phi, E, B, psi, U)
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

!> Calculates the electric and magnetic fields at a specific position
!> in the jorek element `i_elm` at `st`.
subroutine calc_VBpsiU(fields, time, i_elm, st, phi, V, B, psi, U)
use phys_module, only: F0, mode, central_mass, central_density
use constants, only: mu_zero, mass_proton
use mod_coordinate_transforms, only: transform_derivatives_st_to_RZ
! Routine parameters
class(fields_base), intent(in) :: fields
real*8, intent(in)  :: time
integer, intent(in) :: i_elm !< JOREK element index
real*8, intent(in)  :: st(2) !< element-local coordinates
real*8, intent(in)  :: phi !< toroidal angle
real*8, intent(out) :: V(3) !< Parallel velocity in [m/(s)]
real*8, intent(out) :: B(3) !< Magnetic field [T]
real*8, intent(out) :: psi !< psi in JOREK units
real*8, intent(out) :: U !< velocity stream function in m/s

real*8              :: Vpar

! Internal parameters
integer, parameter :: i_var(3) = [var_psi,var_u,var_Vpar]
real*8             :: P(3), P_s(3), P_t(3), P_phi(3), P_time(3) ! Placeholder for evaluating variables and derivatives locally
! Values
real*8             :: R, R_s, R_t, Z, Z_s, Z_t
! Others
real*8             :: inv_st_jac, R_inv
real*8             :: psi_R, psi_Z, U_R, U_Z, U_phi, t_norm

t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

! Interpolate the fields to get psi and U at the current position (and the
! changes u_n - u(n-1))
call fields%interp_PRZ(time, i_elm, i_var, 3, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

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

! Calculate the value of the parallel velocity field
Vpar = P(3)/t_norm

! Calculate the velocity field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
V    = [-U_Z, +U_R, 0.] * R / t_norm
V    = V + Vpar * B

end subroutine calc_VBpsiU

pure subroutine calc_NeTe(fields, time, i_elm, st, phi, n_e, T_e, grad_T_e)
use phys_module, only: central_density
use mod_parameters
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

#if (JOREK_MODEL == 400 || JOREK_MODEL == 502)
! electron temperature
call fields%interp_PRZ(time,i_elm,[5,var_Te],2,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t) 
#else
! electron temperature + ion temperature (assumed equal)
call fields%interp_PRZ(time,i_elm,[5,var_T],2,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
#endif

n_e = max(central_density * P(1) * 1d20,1d16)                           ! plasma density [1/m^3], capped against negative
T_norm = (1.d0/K_BOLTZ/(2.d0*MU_ZERO*central_density*1.d20))
#if (JOREK_MODEL == 400 || JOREK_MODEL == 502)
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

pure subroutine calc_NeTe_imp(fields, time, i_elm, st, phi, n_bg, n_imp, T_e, grad_T_e)
use phys_module, only: central_density, imp_type, central_mass
use constants
use mod_parameters
class(fields_base), intent(in)                    :: fields
integer, intent(in)                               :: i_elm
real*8, intent(in)                                :: time, st(2), phi
real*8, intent(out)                               :: n_bg !< background species density [m^-3]
real*8, intent(out)                               :: n_imp !< Impurity species density [m^-3]
real*8, intent(out)                               :: T_e !< electron temperature [K]
real*8, intent(out), optional, dimension(3)       :: grad_T_e !< gradient of electron temperature [K/m]

real*8, dimension(3) :: P, P_s, P_t, P_phi, P_time
real*8               :: R, R_s, R_t, Z, Z_s, Z_t, xjac, m_i_over_m_imp
real*8 :: T_norm !< temperature normalisation

call fields%interp_PRZ(time,i_elm,&
#if (JOREK_MODEL == 400 || JOREK_MODEL == 502)
      [var_rho,var_Te, var_rhon],& ! electron temperature
#else
      [var_rho,var_T, var_rhon],& ! electron temperature + ion temperature (assumed equal)
#endif
          3,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)

select case ( trim(imp_type) )
  case('D2')
    m_i_over_m_imp = central_mass/2.
  case('Ar')
    m_i_over_m_imp = central_mass/40. ! Argon mass = 40 u and main ion (D) mass = 2 u
  case('Ne')
    m_i_over_m_imp = central_mass/20. ! Neon mass = 20 u and main ion (D) mass = 2 u
  case default
    m_i_over_m_imp = central_mass/2.
end select

n_bg  = P(1) - P(3)
n_imp = P(3) * m_i_over_m_imp
n_bg  = max(central_density * n_bg * 1d20,1d16)        ! plasma density [1/m^3], capped against negative
n_imp = max(central_density * n_imp * 1d20,0.0)        ! plasma density [1/m^3], capped against negative
T_norm = (1.d0/K_BOLTZ/(2.d0*MU_ZERO*central_density*1.d20))
#if (JOREK_MODEL == 400 || JOREK_MODEL == 502)
T_norm = T_norm*2.d0 ! P(1) contains the electron temperature, reverse previous correction
#endif
T_e = max(P(2)*T_norm, 1.d0) ! temperature capped against going negative

if (present(grad_T_e)) then

  xjac = R_s * Z_t - R_t * Z_s
  grad_T_e = T_norm*[(  P_s(2) * Z_t - P_t(2) * Z_s)/ xjac, &
                     (- P_s(2) * R_t + P_t(2) * R_s)/ xjac, &
                     P_phi(2)/R]
end if
end subroutine calc_NeTe_imp

function rot_tmp(x,A,dA) result(rotA)
  implicit none
  real*8, intent(in) :: x(3), A(3), dA(3,3)
  real*8             :: rotA(3)
     rotA(1) = dA(3,2) - dA(2,3) / x(1)
     rotA(2) = dA(1,3) - dA(3,1) - A(3) / x(1)
     rotA(3) = dA(2,1) - dA(1,2) 
  return
end  

subroutine calc_RK4_analytic(fields, R, Z, phi, A_out, dA_out, B_out, dB_out, B_norm, dB_norm, bn, dBn, E)
  use phys_module, only: mode, central_mass, central_density
  use constants, only: mu_zero, mass_proton
  class(fields_base), intent(in) :: fields
  ! Routine parameters
  real*8, intent(in)  :: R,Z, phi      !< position in  [m,m,rad]
  real*8, intent(out) :: E(3)          !< Electric field [V/m]
  real*8, intent(out) :: A_out(3)      !< vector potential [T.m]
  real*8, intent(out) :: dA_out(3,3)   !< derivatives of vector potential [T]
  real*8, intent(out) :: B_out(3)      !< Magnetic field [T]
  real*8, intent(out) :: dB_out(3,3)   !< derivatives of magnetic field [T/m]
  real*8, intent(out) :: B_norm(3)     !< normalised magnetic field vector
  real*8, intent(out) :: dB_norm(3,3)  !< derivatives of normalised magnetic field vector [1/m]
  real*8, intent(out) :: Bn            !< Magnetic field amplitude [T]
  real*8, intent(out) :: dBn(3)        !< derivatives of magnetic field amplitude [T/m]

  real*8 :: AR, AR_R, AR_Z, AR_phi
  real*8 :: AZ, AZ_R, AZ_Z, AZ_phi
  real*8 :: Aphi, Aphi_R, Aphi_Z, Aphi_phi
  real*8 :: BR, BR_R, BR_Z, BR_phi
  real*8 :: BZ, BZ_R, BZ_Z, BZ_phi
  real*8 :: Bphi, BPhi_R, Bphi_Z, Bphi_phi
  real*8 :: BBR, BBZ, BBphi
  real*8 :: R0, B0, F0, q, S, dS_R, dS_Z, dS_phi, rr, rr_R, rr_Z

  R0 = 1.d0
  B0 = 1.d0
  F0 = 1.d0
  q  = 2

  AR     =   F0 * Z / (2.d0 * R)
  AR_R   = - F0 * Z / (2.d0 * R**2)
  AR_Z   = + F0     / (2.d0 * R)
  AR_phi = 0.d0

  AZ     = - log(R/R0) * F0 /2.d0
  AZ_R   = - 1.d0/R    * F0 /2.d0
  AZ_Z   = 0.d0
  AZ_phi = 0.d0

  rr   = sqrt((R-R0)**2 + Z**2)    
  rr_R = 1.d0 / (2.d0 * rr) * 2.d0*(R-R0)
  rr_Z = 1.d0 / (2.d0 * rr) * 2.d0*Z

  Aphi   = - B0 * rr**2 / (2.d0 * q * R)
  Aphi_R = - B0 * 2.d0 * rr * rr_R / (2.d0 * q * R) + B0 * rr**2 / (2.d0 * q * R**2)
  Aphi_Z = - B0 * 2.d0 * rr * rr_Z / (2.d0 * q * R)
  Aphi_phi = 0.d0

  BR     = - B0 * Z / (q * R)
  BR_R   = + B0 * Z / (q * R**2)
  BR_Z   = - B0     / (q * R)
  BR_phi = 0.d0

  BZ     = B0*(R-R0) / (q * R)
  BZ_R   = B0 / (q * R) -  B0*(R-R0) / (q * R**2)
  BZ_Z   = 0.d0
  BZ_phi = 0.d0 

  Bphi     = - F0 / R
  Bphi_R   = F0 / R**2
  Bphi_Z   = 0.d0
  Bphi_phi = 0.d0

  BBR   = Aphi_Z - AZ_phi / R
  BBZ   = AR_phi - Aphi/R - Aphi_R
  BBphi = AZ_R   - AR_Z 

  S      = sqrt((R - R0)**2 + Z**2 + q*2 * R0**2)
  dS_R   = 1.d0/(2.d0*S) * 2.d0*(R - R0)
  dS_Z   = 1.d0/(2.d0*S) * 2.d0*z
  dS_phi = 0.d0

  Bn     = B0 * S / (q * R)
  dBn(1) = B0 * dS_R   / (q * R) - B0 * S /(q * R**2)
  dBn(2) = B0 * dS_Z   / (q * R) 
  dBn(3) = B0 * dS_phi / (q * R) 
  
  A_out(1) = AR;    A_out(2) = AZ;    A_out(3) = Aphi
  B_out(1) = BR;    B_out(2) = BZ;    B_out(3) = Bphi

  dA_out(1,1) = AR_R;    dA_out(1,2) = AR_Z;    dA_out(1,3) = AR_phi
  dA_out(2,1) = AZ_R;    dA_out(2,2) = AZ_Z;    dA_out(2,3) = AZ_phi
  dA_out(3,1) = Aphi_R;  dA_out(3,2) = Aphi_Z;  dA_out(3,3) = Aphi_phi

  dB_out(1,1) = BR_R;    dB_out(1,2) = BR_Z;    dB_out(1,3) = BR_phi
  dB_out(2,1) = BZ_R;    dB_out(2,2) = BZ_Z;    dB_out(2,3) = BZ_phi
  dB_out(3,1) = Bphi_R;  dB_out(3,2) = Bphi_Z;  dB_out(3,3) = Bphi_phi

  B_norm = B_out / Bn

  dB_norm(1,:) = dB_out(1,:) / Bn - B_out(1) / Bn**2 * dBn(:)
  dB_norm(2,:) = dB_out(2,:) / Bn - B_out(2) / Bn**2 * dBn(:)
  dB_norm(3,:) = dB_out(3,:) / Bn - B_out(3) / Bn**2 * dBn(:)
  
  E = 0.d0

return
end

subroutine calc_RK4(fields, time, i_elm, st, phi, A, dA, B, dB, Bnorm, dBnorm, bn, dBn, E)
use phys_module, only: F0, mode, central_mass, central_density
use constants, only: mu_zero, mass_proton
! Routine parameters
class(fields_base), intent(in) :: fields
real*8, intent(in)  :: time
integer, intent(in) :: i_elm       !< JOREK element index
real*8, intent(in)  :: st(2)       !< element-local coordinates
real*8, intent(in)  :: phi         !< toroidal angle
real*8, intent(out) :: E(3)        !< Electric field [V/m]
real*8, intent(out) :: A(3)        !< vector potential [T.m]
real*8, intent(out) :: dA(3,3)     !< derivatives of vector potential [T]
real*8, intent(out) :: B(3)        !< Magnetic field [T]
real*8, intent(out) :: dB(3,3)     !< derivatives of magnetic field [T/m]
real*8, intent(out) :: Bnorm(3)    !< normalised magnetic field vector
real*8, intent(out) :: dBnorm(3,3) !< derivatives of normalised magnetic field vector [1/m]
real*8, intent(out) :: Bn          !< Magnetic field amplitude [T]
real*8, intent(out) :: dBn(3)      !< derivatives of magnetic field amplitude [T/m]

! Internal parameters
integer, parameter :: i_var(2) = [1,2]
real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_ss(2), P_st(2), P_tt(2), P_sphi(2), P_tphi(2)
real*8             :: P_time(2), P_stime(2), P_ttime(2), bn2
real*8             :: x(3), R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
! Others
real*8             :: inv_st_jac, R_inv, RZjac, RZjac_R, RZjac_Z
real*8             :: psi, psi_R, psi_Z, psi_RR, psi_ZZ, psi_RZ, psi_Rphi, psi_Zphi
real*8             :: U, U_R, U_Z, U_phi, t_norm

t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

call fields%interp_PRZ_2(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, &
                       P_time, P_ss, P_st, P_tt, P_sphi, P_tphi, P_stime, P_ttime,   &
                       R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_ss, Z_st, Z_tt)

R_inv = 1.d0/R
inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)

! Update psi and U
psi = P(1)
U   = P(2)/t_norm

! Calculate the derivatives to R and Z
psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
U_R      = (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
U_Z      = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
U_phi    = P_phi(2)

psi_Rphi = (  P_sphi(1) * Z_t - P_tphi(1) * Z_s ) * inv_st_jac
psi_Zphi = (- P_sphi(1) * R_t + P_tphi(1) * R_s ) * inv_st_jac

RZjac    = R_s*Z_t - R_t*Z_s

RZjac_R  = (R_ss*Z_t**2 - Z_ss*R_t*Z_t - 2.d0*R_st*Z_s*Z_t   &
         + Z_st*(R_s*Z_t + R_t*Z_s) + R_tt*Z_s**2 - Z_tt*R_s*Z_s) / RZjac

RZjac_Z  = (Z_tt*R_s**2 - R_tt*Z_s*R_s - 2.d0*Z_st*R_t*R_s   &
         + R_st*(Z_t*R_s + Z_s*R_t) + Z_ss*R_t**2 - R_ss*Z_t*R_t) / RZjac

psi_RR = (P_ss(1) * Z_t**2 - 2.d0*P_st(1) * Z_s*Z_t + P_tt(1) * Z_s**2               &
       + P_s(1) * (Z_st*Z_t - Z_tt*Z_s) + P_t(1) * (Z_st*Z_s - Z_ss*Z_t)) / RZjac**2 &
       - RZjac_R * (P_s(1) * Z_t - P_t(1) * Z_s) / RZjac**2

psi_ZZ = (P_ss(1) * R_t**2 - 2.d0*P_st(1) * R_s*R_t + P_tt(1) * R_s**2                &
       + P_s(1) * (R_st*R_t - R_tt*R_s ) + P_t(1) * (R_st*R_s - R_ss*R_t)) / RZjac**2 &
       - RZjac_Z * (- P_s(1) * R_t + P_t(1) * R_s) / RZjac**2

psi_RZ = (- P_ss(1) * Z_t*R_t - P_tt(1) * R_s*Z_s + P_st(1) * (Z_s*R_t + Z_t*R_s)       &
       - P_s(1) * (R_st*Z_t - R_tt*Z_s) - P_t(1) * (R_st*Z_s - R_ss*Z_t) )  / RZjac**2  &
       - RZjac_R * (- P_s(1) * R_t + P_t(1) * R_s)   / RZjac**2

x(1) = R
x(2) = Z
x(3) = phi

A = (/ - F0 * Z / (2.d0 * R),  + log(R) * F0 /2.d0, psi / R /)

dA(1,1) = + F0 * Z / (2.d0 * R**2)
dA(1,2) = - F0     / (2.d0 * R)
dA(1,3) = 0.d0
dA(2,1) = + 1.d0/R * F0 /2.d0
dA(2,2) = 0.d0
dA(2,3) = 0.d0
dA(3,1) = psi_R / R - psi / R**2
dA(3,2) = psi_Z / R
dA(3,3) = P_phi(1) / R

! Set dpsi/dt to 0 if flag is true
if(fields%flag_zero_dpsidt) P_time(1) = 0.d0

! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
B     = [+psi_Z, -psi_R, F0] * R_inv

dB(1,1) =   psi_RZ   * R_inv - psi_Z * R_inv**2
dB(1,2) =   psi_ZZ   * R_inv
dB(1,3) =   psi_Zphi * R_inv
dB(2,1) = - psi_RR   * R_inv + psi_R * R_inv**2
dB(2,2) = - psi_RZ   * R_inv
dB(2,3) = - psi_Rphi * R_inv
dB(3,1) = - F0 / R**2 ! additional terms for toroidal geometry?
dB(3,2) =  0.d0
dB(3,3) =  0.d0

Bn    = norm2(B)
Bn2   = Bn**2
Bnorm = B / Bn

Bn = sqrt(psi_R**2 + psi_Z**2 + F0**2) / R

dBn(1) = 0.5d0 /(R*Bn) * (2.d0*psi_R * psi_RR   + 2.d0*psi_Z * psi_RZ) - Bn / R
dBn(2) = 0.5d0 /(R*Bn) * (2.d0*psi_R * psi_RZ   + 2.d0*psi_Z * psi_ZZ)
dBn(3) = 0.5d0 /(R*Bn) * (2.d0*psi_R * psi_Rphi + 2.d0*psi_Z * psi_Zphi)

dBnorm(1,1) = dB(1,1) / Bn - B(1) * dBn(1) / Bn2
dBnorm(1,2) = dB(1,2) / Bn - B(1) * dBn(2) / Bn2
dBnorm(1,3) = dB(1,3) / Bn - B(1) * dBn(3) / Bn2
dBnorm(2,1) = dB(2,1) / Bn - B(2) * dBn(1) / Bn2
dBnorm(2,2) = dB(2,2) / Bn - B(2) * dBn(2) / Bn2
dBnorm(2,3) = dB(2,3) / Bn - B(2) * dBn(3) / Bn2
dBnorm(3,1) = dB(3,1) / Bn - B(3) * dBn(1) / Bn2
dBnorm(3,2) = dB(3,2) / Bn - B(3) * dBn(2) / Bn2
dBnorm(3,3) = dB(3,3) / Bn - B(3) * dBn(3) / Bn2

dBnorm(3,:) = dBnorm(3,:) * R_inv

! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
! See http://jorek.eu/wiki/doku.php?id=u_phi
E     = [-F0*U_R, -F0*U_Z, -F0*U_phi*R_inv]/t_norm
E(3)  = E(3) - R_inv*P_time(1) ! because this is not normalized with t_norm

end subroutine calc_RK4

subroutine calc_Qin_analytic(fields, R, Z, phi, A_out, dA_out, B_out, dB_out, B_norm, dB_norm, bn, dBn, E)
  use phys_module, only: mode, central_mass, central_density
  use constants, only: mu_zero, mass_proton
  class(fields_base), intent(in) :: fields
  ! Routine parameters
  real*8, intent(in)  :: R,Z, phi      !< position in  [m,m,rad]
  real*8, intent(out) :: E(3)          !< Electric field [V/m]
  real*8, intent(out) :: A_out(3)      !< vector potential [T.m]
  real*8, intent(out) :: dA_out(3,3)   !< derivatives of vector potential [T]
  real*8, intent(out) :: B_out(3)      !< Magnetic field [T]
  real*8, intent(out) :: dB_out(3,3)   !< derivatives of magnetic field [T/m]
  real*8, intent(out) :: B_norm(3)     !< normalised magnetic field vector
  real*8, intent(out) :: dB_norm(3,3)  !< derivatives of normalised magnetic field vector [1/m]
  real*8, intent(out) :: Bn            !< Magnetic field amplitude [T]
  real*8, intent(out) :: dBn(3)        !< derivatives of magnetic field amplitude [T/m]

  real*8 :: AR, AR_R, AR_Z, AR_phi
  real*8 :: AZ, AZ_R, AZ_Z, AZ_phi
  real*8 :: Aphi, Aphi_R, Aphi_Z, Aphi_phi
  real*8 :: BR, BR_R, BR_Z, BR_phi
  real*8 :: BZ, BZ_R, BZ_Z, BZ_phi
  real*8 :: Bphi, BPhi_R, Bphi_Z, Bphi_phi
  real*8 :: BBR, BBZ, BBphi
  real*8 :: R0, B0, F0, q, S, dS_R, dS_Z, dS_phi, rr, rr_R, rr_Z

  R0 = 1.d0
  B0 = 1.d0
  F0 = 1.d0
  q  = 2

  E = 0.d0

  AR     =   F0 * Z / (2.d0 * R)
  AR_R   = - F0 * Z / (2.d0 * R**2)
  AR_Z   = + F0     / (2.d0 * R)
  AR_phi = 0.d0

  AZ     = - log(R/R0) * F0 /2.d0
  AZ_R   = - 1.d0/R    * F0 /2.d0
  AZ_Z   = 0.d0
  AZ_phi = 0.d0

  rr   = sqrt((R-R0)**2 + Z**2)    
  rr_R = 1.d0 / (2.d0 * rr) * 2.d0*(R-R0)
  rr_Z = 1.d0 / (2.d0 * rr) * 2.d0*Z

  Aphi   = - B0 * rr**2 / (2.d0 * q * R)
  Aphi_R = - B0 * 2.d0 * rr * rr_R / (2.d0 * q * R) + B0 * rr**2 / (2.d0 * q * R**2)
  Aphi_Z = - B0 * 2.d0 * rr * rr_Z / (2.d0 * q * R)
  Aphi_phi = 0.d0

  BR     = - B0 * Z / (q * R)
  BR_R   = + B0 * Z / (q * R**2)
  BR_Z   = - B0     / (q * R)
  BR_phi = 0.d0

  BZ     = B0*(R-R0) / (q * R)
  BZ_R   = B0 / (q * R) -  B0*(R-R0) / (q * R**2)
  BZ_Z   = 0.d0
  BZ_phi = 0.d0 

  Bphi     = - F0 / R
  Bphi_R   = F0 / R**2
  Bphi_Z   = 0.d0
  Bphi_phi = 0.d0

  BBR   = Aphi_Z - AZ_phi / R
  BBZ   = AR_phi - Aphi/R - Aphi_R
  BBphi = AZ_R   - AR_Z 

  S      = sqrt((R - R0)**2 + Z**2 + q*2 * R0**2)
  dS_R   = 1.d0/(2.d0*S) * 2.d0*(R - R0)
  dS_Z   = 1.d0/(2.d0*S) * 2.d0*z
  dS_phi = 0.d0

  Bn     = B0 * S / (q * R)
  dBn(1) = B0 * dS_R   / (q * R) - B0 * S /(q * R**2)
  dBn(2) = B0 * dS_Z   / (q * R) 
  dBn(3) = B0 * dS_phi / (q * R) 
  
  A_out(1) = AR;    A_out(2) = AZ;    A_out(3) = Aphi
  B_out(1) = BR;    B_out(2) = BZ;    B_out(3) = Bphi

  dA_out(1,1) = AR_R;    dA_out(1,2) = AR_Z;    dA_out(1,3) = AR_phi
  dA_out(2,1) = AZ_R;    dA_out(2,2) = AZ_Z;    dA_out(2,3) = AZ_phi
  dA_out(3,1) = Aphi_R;  dA_out(3,2) = Aphi_Z;  dA_out(3,3) = Aphi_phi

  dB_out(1,1) = BR_R;    dB_out(1,2) = BR_Z;    dB_out(1,3) = BR_phi
  dB_out(2,1) = BZ_R;    dB_out(2,2) = BZ_Z;    dB_out(2,3) = BZ_phi
  dB_out(3,1) = Bphi_R;  dB_out(3,2) = Bphi_Z;  dB_out(3,3) = Bphi_phi

  B_norm = B_out / Bn

  dB_norm(1,:) = dB_out(1,:) / Bn - B_out(1) / Bn**2 * dBn(:)
  dB_norm(2,:) = dB_out(2,:) / Bn - B_out(2) / Bn**2 * dBn(:)
  dB_norm(3,:) = dB_out(3,:) / Bn - B_out(3) / Bn**2 * dBn(:)

  !----------------- convert to covariant toroidal component
  dA_out(3,1)  = R * dA_out(3,1) + A_out(3)
  dA_out(3,2)  = R * dA_out(3,2)
  dA_out(3,3)  = R * dA_out(3,3)

  dB_norm(3,1) = R * dB_norm(3,1) + B_norm(3)
  dB_norm(3,2) = R * dB_norm(3,2)
  dB_norm(3,3) = R * dB_norm(3,3)

  dB_out(3,1) = R * dB_out(3,1) + B_out(3)
  dB_out(3,2) = R * dB_out(3,2)
  dB_out(3,3) = R * dB_out(3,3)

  A_out(3)     = R * A_out(3)
  B_norm(3)    = R * B_norm(3)
  B_out(3)     = R * B_out(3)

  dBn(3) = R * dBn(3)
  E(3)  = R * E(3)

return
end

subroutine calc_Qin(fields, time, i_elm, st, phi, A, dA, B, dB, Bnorm, dBnorm, bn, dBn, E)
  use phys_module, only: F0, mode, central_mass, central_density
  use constants, only: mu_zero, mass_proton
  ! Routine parameters
  class(fields_base), intent(in) :: fields
  real*8, intent(in)  :: time
  integer, intent(in) :: i_elm       !< JOREK element index
  real*8, intent(in)  :: st(2)       !< element-local coordinates
  real*8, intent(in)  :: phi         !< toroidal angle
  real*8, intent(out) :: E(3)        !< Electric field [V/m]
  real*8, intent(out) :: A(3)        !< vector potential [T.m]
  real*8, intent(out) :: dA(3,3)     !< derivatives of vector potential [T]
  real*8, intent(out) :: B(3)        !< Magnetic field [T]
  real*8, intent(out) :: dB(3,3)     !< derivatives of magnetic field [T/m]
  real*8, intent(out) :: Bnorm(3)    !< normalised magnetic field vector
  real*8, intent(out) :: dBnorm(3,3) !< derivatives of normalised magnetic field vector [1/m]
  real*8, intent(out) :: Bn          !< Magnetic field amplitude [T]
  real*8, intent(out) :: dBn(3)      !< derivatives of magnetic field amplitude [T/m]
  
  ! Internal parameters
  integer, parameter :: i_var(2) = [1,2]
  real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_ss(2), P_st(2), P_tt(2), P_sphi(2), P_tphi(2)
  real*8             :: P_time(2), P_stime(2), P_ttime(2), bn2
  real*8             :: x(3), R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
  ! Others
  real*8             :: inv_st_jac, R_inv, RZjac, RZjac_R, RZjac_Z
  real*8             :: psi, psi_R, psi_Z, psi_RR, psi_ZZ, psi_RZ, psi_Rphi, psi_Zphi
  real*8             :: U, U_R, U_Z, U_phi, t_norm
  
  t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds
  
  call fields%interp_PRZ_2(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, &
                         P_time, P_ss, P_st, P_tt, P_sphi, P_tphi, P_stime, P_ttime,   &
                         R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_ss, Z_st, Z_tt)

  R_inv = 1.d0/R
  inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
    
  ! Update psi and U
  psi = P(1)
  U   = P(2)/t_norm
  
  ! Calculate the derivatives to R and Z
  psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
  psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
  U_R      = (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
  U_Z      = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
  U_phi    = P_phi(2)
  
  psi_Rphi = (  P_sphi(1) * Z_t - P_tphi(1) * Z_s ) * inv_st_jac
  psi_Zphi = (- P_sphi(1) * R_t + P_tphi(1) * R_s ) * inv_st_jac
  
  RZjac    = R_s*Z_t - R_t*Z_s
  
  RZjac_R  = (R_ss*Z_t**2 - Z_ss*R_t*Z_t - 2.d0*R_st*Z_s*Z_t   &
           + Z_st*(R_s*Z_t + R_t*Z_s) + R_tt*Z_s**2 - Z_tt*R_s*Z_s) / RZjac
  
  RZjac_Z  = (Z_tt*R_s**2 - R_tt*Z_s*R_s - 2.d0*Z_st*R_t*R_s   &
           + R_st*(Z_t*R_s + Z_s*R_t) + Z_ss*R_t**2 - R_ss*Z_t*R_t) / RZjac
  
  psi_RR = (P_ss(1) * Z_t**2 - 2.d0*P_st(1) * Z_s*Z_t + P_tt(1) * Z_s**2               &
         + P_s(1) * (Z_st*Z_t - Z_tt*Z_s) + P_t(1) * (Z_st*Z_s - Z_ss*Z_t)) / RZjac**2 &
         - RZjac_R * (P_s(1) * Z_t - P_t(1) * Z_s) / RZjac**2
  
  psi_ZZ = (P_ss(1) * R_t**2 - 2.d0*P_st(1) * R_s*R_t + P_tt(1) * R_s**2                &
         + P_s(1) * (R_st*R_t - R_tt*R_s ) + P_t(1) * (R_st*R_s - R_ss*R_t)) / RZjac**2 &
         - RZjac_Z * (- P_s(1) * R_t + P_t(1) * R_s) / RZjac**2
  
  psi_RZ = (- P_ss(1) * Z_t*R_t - P_tt(1) * R_s*Z_s + P_st(1) * (Z_s*R_t + Z_t*R_s)       &
         - P_s(1) * (R_st*Z_t - R_tt*Z_s) - P_t(1) * (R_st*Z_s - R_ss*Z_t) )  / RZjac**2  &
         - RZjac_R * (- P_s(1) * R_t + P_t(1) * R_s)   / RZjac**2
  
  x(1) = R
  x(2) = Z
  x(3) = phi
  
  A = (/ - F0 * Z / (2.d0 * R),  + log(R) * F0 /2.d0, psi / R /)
  
  dA(1,1) = + F0 * Z / (2.d0 * R**2)
  dA(1,2) = - F0     / (2.d0 * R)
  dA(1,3) = 0.d0
  dA(2,1) = + 1.d0/R * F0 /2.d0
  dA(2,2) = 0.d0
  dA(2,3) = 0.d0
  dA(3,1) = psi_R / R - psi / R**2
  dA(3,2) = psi_Z / R
  dA(3,3) = P_phi(1) / R
  
  ! Set dpsi/dt to 0 if flag is true
  if(fields%flag_zero_dpsidt) P_time(1) = 0.d0
  
  ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
  B     = [+psi_Z, -psi_R, F0] * R_inv
  
  dB(1,1) =   psi_RZ   * R_inv - psi_Z * R_inv**2
  dB(1,2) =   psi_ZZ   * R_inv
  dB(1,3) =   psi_Zphi * R_inv
  dB(2,1) = - psi_RR   * R_inv + psi_R * R_inv**2
  dB(2,2) = - psi_RZ   * R_inv
  dB(2,3) = - psi_Rphi * R_inv
  dB(3,1) = - F0 / R**2 ! additional terms for toroidal geometry?
  dB(3,2) =  0.d0
  dB(3,3) =  0.d0
  
  Bn    = norm2(B)
  Bn2   = Bn**2
  Bnorm = B / Bn
  
  Bn = sqrt(psi_R**2 + psi_Z**2 + F0**2) / R
  
  dBn(1) = 0.5d0 /(R*Bn) * (2.d0*psi_R * psi_RR   + 2.d0*psi_Z * psi_RZ) - Bn / R
  dBn(2) = 0.5d0 /(R*Bn) * (2.d0*psi_R * psi_RZ   + 2.d0*psi_Z * psi_ZZ)
  dBn(3) = 0.5d0 /(R*Bn) * (2.d0*psi_R * psi_Rphi + 2.d0*psi_Z * psi_Zphi)
  
  dBnorm(1,:) = dB(1,:) / Bn - B(1) * dBn(:) / Bn2
  dBnorm(2,:) = dB(2,:) / Bn - B(2) * dBn(:) / Bn2
  dBnorm(3,:) = dB(3,:) / Bn - B(3) * dBn(:) / Bn2
    
  ! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
  ! See http://jorek.eu/wiki/doku.php?id=u_phi
  E     = [-F0*U_R, -F0*U_Z, -F0*U_phi*R_inv]/t_norm
  E(3)  = E(3) - R_inv*P_time(1) ! because this is not normalized with t_norm

  !----------------- convert to covariant toroidal component
  dA(3,1)  = R * dA(3,1) + A(3)
  dA(3,2)  = R * dA(3,2)
  dA(3,3)  = R * dA(3,3) 

  dBnorm(3,1) = R * dBnorm(3,1) + Bnorm(3)
  dBnorm(3,2) = R * dBnorm(3,2)
  dBnorm(3,3) = R * dBnorm(3,3)

  dB(3,1) = R * dB(3,1) + B(3)
  dB(3,2) = R * dB(3,2)
  dB(3,3) = R * dB(3,3)

  A(3)     = R * A(3)
  Bnorm(3) = R * Bnorm(3)
  B(3)     = R * B(3)
  dBn(3)   = R * dBn(3)
  E(3)     = R * E(3)
  
end subroutine calc_Qin
  

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
