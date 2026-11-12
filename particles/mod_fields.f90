!> Module containing base type for field interpolations, with interfaces
!> to implement
module mod_fields
  use data_structure
  use corr_neg
  
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
    procedure(interp_PRZ), deferred, public    :: interp_PRZ
    procedure(interp_PRZ_2), deferred, public  :: interp_PRZ_2
    procedure(interp_PRZP_1), deferred, public :: interp_PRZP_1
    procedure, public :: calc_NeTe
    procedure, public :: calc_vpar
    procedure, public :: calc_NeTeTi
    procedure, public :: calc_NjTj
    procedure, public :: calc_EBpsiU
    procedure, public :: calc_vvector
    procedure, public :: calc_F_profile
    procedure, public :: calc_gyro_average_E
    procedure, public :: calc_Qin, calc_Qin_analytic, check_consistency_Qin
    procedure, public :: calc_rk4, calc_RK4_analytic, check_consistency_RK4
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
    pure subroutine interp_PRZ_2(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, &
                               P_time, P_ss, P_st, P_tt, P_sphi, P_tphi, P_stime, P_ttime, &
                               R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_ss, Z_st, Z_tt)
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
    !> Interpolate a variable at s, t, phi in i_elm, returning first
    !> derivatives of the variable and of space
    pure subroutine interp_PRZP_1(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, R_phi, Z, Z_s, Z_t, Z_phi)
      import fields_base
      class(fields_base),  intent(in)  :: this
      real*8,                   intent(in)  :: time !< Time at which to calculate this variable
      integer,                  intent(in)  :: i_elm
      integer,                  intent(in)  :: n_v, i_v(n_v)
      real*8,                   intent(in)  :: s, t, phi
      real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v), P_time(n_v)
      real*8,                   intent(out) :: R, R_s, R_t, R_phi, Z, Z_s, Z_t, Z_phi
    end subroutine interp_PRZP_1
  end interface

contains
!> Calculates the electric and magnetic fields at a specific position
!> in the jorek element `i_elm` at `st`.
subroutine calc_EBpsiU(fields, time, i_elm, st, phi, E, B, psi, U)
  use phys_module, only: F0, mode, central_mass, central_density
  use constants, only: mu_zero, atomic_mass_unit
  use mod_coordinate_transforms, only: transform_derivatives_st_to_RZ
  use mod_chi
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
#ifdef fullmhd
  integer, parameter :: i_var(3) = [1,2,3]
  real*8             :: P(3), P_s(3), P_t(3), P_phi(3), P_time(3) !Placeholders, differ in full MHD
  real*8             :: A3, AR, AZ, A3_R, A3_Z, A3_t, AR_Z, AR_p, AR_t, AZ_R, AZ_P,AZ_t, Fprof
#else
  integer, parameter :: i_var(2) = [1,2]
  real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_time(2) ! Placeholder for evaluating variables and derivatives locally
#endif
  ! Values
  real*8             :: R, R_s, R_t, R_phi, Z, Z_s, Z_t, Z_phi
  ! Others
  real*8             :: inv_st_jac, R_inv
  real*8             :: psi_R, psi_Z, psi_phi, U_R, U_Z, U_phi, t_norm
  real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

  t_norm  = sqrt(mu_zero * ATOMIC_MASS_UNIT * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  ! Interpolate the fields to get psi and U at the current position (and the
  ! changes u_n - u(n-1))


#ifdef fullmhd
  call fields%calc_F_profile(i_elm,st(1),st(2),phi,Fprof)
  !In full MHD equations are easier,
  ! B = F/R e_\phi  + curl A
  ! E=\partial_t A
  !
  !Interpolating A^3=psi=1, AR=2, AZ=3, including time derivatives
  call fields%interp_PRZ(time, i_elm, i_var,3, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

  R_inv = 1.d0/R
  inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
  A3=P(1)
  AR=P(2)
  AZ=P(3)
  !Derivatives of A3
  A3_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
  A3_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
  A3_t    = P_time(1)
  !Derivatives of AR
  AR_Z    = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
  AR_p    = P_phi(2)
  AR_t    = P_time(2)
  !Derivatives of AZ
  AZ_R    = (  P_s(3) * Z_t - P_t(3) * Z_s ) * inv_st_jac
  AZ_p    = P_phi(3)
  AZ_t    = P_time(3)

  B=[(A3_Z-AZ_p)*R_inv, (AR_p-A3_R)*R_inv, AZ_R-AR_Z + Fprof*R_inv]
  E=[-AR_t, -AZ_t, -R_inv*A3_t]
#else
#if STELLARATOR_MODEL
  call fields%interp_PRZP_1(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, R_phi, Z, Z_s, Z_t, Z_phi)
#else
  R_phi = 0.0; Z_phi = 0.0
  call fields%interp_PRZ(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
#endif
  ! Calculate the derivatives to R and Z

  R_inv = 1.d0/R
  inv_st_jac = 1.d0/jac(R_s,R_t,Z_s,Z_t)
  psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
  psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
  psi_phi  = P_phi(1) - R_phi*psi_R - Z_phi*psi_Z
  U_R      = (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
  U_Z      = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
  U_phi    = P_phi(2) - R_phi*U_R - Z_phi*U_Z

  ! Update psi and U
  psi = P(1)
  U   = P(2)/t_norm

  ! Set dpsi/dt to 0 if flag is true
  if(fields%flag_zero_dpsidt) P_time(1) = 0.d0

#if STELLARATOR_MODEL
  chi = get_chi(R,Z,phi,fields%node_list,fields%element_list,i_elm,st(1),st(2))

  B = (/ chi(1,0,0)      + (psi_Z*chi(0,0,1) - psi_phi*chi(0,1,0))/(F0*R), &
         chi(0,1,0)      - (psi_R*chi(0,0,1) - psi_phi*chi(1,0,0))/(F0*R), &
         chi(0,0,1)/R    + (psi_R*chi(0,1,0) - psi_Z * chi(1,0,0))/F0         /)
  
  E     = [-U_R, -U_Z, -U_phi*R_inv]/t_norm
  E     = E - P_time(1)/F0*(/ chi(1,0,0), chi(0,1,0), chi(0,0,1)*R_inv /) 
#else
  ! Calculate the magnetic field (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
  B     = [+psi_Z, -psi_R, F0] * R_inv

  ! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
  ! See http://jorek.eu/wiki/doku.php?id=u_phi
  E     = [-F0*U_R, -F0*U_Z, -F0*U_phi*R_inv]/t_norm
  E(3)  = E(3) - R_inv*P_time(1) ! because this is not normalized with t_norm
#endif

#endif

end subroutine calc_EBpsiU

pure subroutine calc_F_profile(fields,i_elm,s,t,phi,Fprof)
  use data_structure
  use phys_module, only : mode, F0
  use mod_basisfunctions
  class(fields_base),         intent(in)     :: fields
  integer,                    intent(in)     :: i_elm
  real*8,                     intent(in)     :: s,t, phi
  real*8,                     intent(out)    :: Fprof
  !Internal variables
  integer           :: i,j,i_tor, iv, i_harm
  real*8            :: Fprof_temp
  real*8            :: H(n_vertex_max,n_degrees), H_s(n_vertex_max,n_degrees),H_t(n_vertex_max,n_degrees),ss
#ifdef fullmhd
  Fprof_temp=0.d0
  call basisfunctions3(s,t,H,H_s,H_t)
  do i = 1, n_vertex_max
    iv=fields%element_list%element(i_elm)%vertex(i)
    do j=1, n_degrees
      ss=fields%element_list%element(i_elm)%size(i,j)
      Fprof_temp = Fprof_temp +fields%node_list%node(iv)%Fprof_eq(j)*ss*H(i,j)
    enddo!order
  enddo  !number of vertices
  Fprof=Fprof_temp
#else
  Fprof=F0 !In reduced MHD, there is no F profile.
#endif

end subroutine calc_F_profile

!> Wrapper around the general Te/Ti routine to avoid code duplication.
!> We use keyword arguments to map inputs correctly and skip Ti-related outputs.
subroutine calc_NeTe(fields, time, i_elm, st, phi, n_e, T_e, n_e_raw, T_e_raw, grad_T_e)
  class(fields_base), intent(in)                    :: fields
  integer, intent(in)                               :: i_elm
  real*8, intent(in)                                :: time, st(2), phi
  real*8, intent(out)                               :: n_e 
  real*8, intent(out)                               :: T_e 
  real*8, intent(out), optional                     :: n_e_raw 
  real*8, intent(out), optional                     :: T_e_raw 
  real*8, intent(out), optional, dimension(3)       :: grad_T_e 

  call fields%calc_NeTeTi(time, i_elm, st, phi, &
                          n_e      = n_e,       &
                          T_e      = T_e,       &
                          n_e_raw  = n_e_raw,   &
                          T_e_raw  = T_e_raw,   &
                          grad_T_e = grad_T_e)

end subroutine calc_NeTe


!> Calculates electron density and temperatures (Te and optional Ti) at a specific position
!> (element i_elm, local coords st, toroidal angle phi) ands time.
!> Returns corrected SI values and, if requested, raw values and stemperature gradients,
!> supports both one- and two-temperature models.
subroutine calc_NeTeTi(fields,time,i_elm,st,phi,                   &
                            n_e,T_e,T_i,                                &
                            n_e_raw,T_e_raw,T_i_raw,                    &
                            grad_raw_n_e,grad_T_e,grad_T_i)
  use phys_module, only: central_density
  use corr_neg, only: corr_neg_dens, corr_neg_temp
  use constants
  class(fields_base), intent(in)                    :: fields
  integer, intent(in)                               :: i_elm
  real*8, intent(in)                                :: time, st(2), phi
  real*8, intent(out)                               :: n_e                  !< corrected ne [m^-3]
  real*8, intent(out)                               :: T_e                  !< corrected Te [K]
  real*8, intent(out),  optional                    :: T_i                  !< corrected Ti [K]
  real*8, intent(out),  optional                    :: n_e_raw              !< raw ne [m^-3]
  real*8, intent(out),  optional                    :: T_e_raw, T_i_raw     !< raw Ti/Te [K]
  real*8, intent(out),  optional, dimension(3)      :: grad_raw_n_e         !< grad raw ne [m^-4]
  real*8, intent(out),  optional, dimension(3)      :: grad_T_e, grad_T_i   !< grad(corrected Ti/Te) [K/m]
  
#ifdef WITH_TiTe
  real*8, dimension(3) :: P, P_s, P_t, P_phi, P_time
#else
  real*8, dimension(2) :: P, P_s, P_t, P_phi, P_time
#endif
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t, xjac
  real*8               :: inv_xjac, inv_R
  real*8               :: T_norm, n_norm, Te_jor, Ti_jor
  real*8               :: tmp_g(3)
  integer              :: ii_Ti, ii_Te
  logical              :: need_Ti, need_grad
  real*8, parameter    :: EPS     = 1.d-12

  ! normalizations
  T_norm = 1.d0 / (K_BOLTZ * MU_ZERO * central_density * 1.d20)
  n_norm = central_density * 1.d20

  ! what do we actually need?
  need_Ti   = present(T_i) .or. present(T_i_raw) .or. present(grad_T_i)
  need_grad = present(grad_T_i) .or. present(grad_T_e) .or. present(grad_raw_n_e)

  ! interpolate fields
#ifdef WITH_TiTe
  call fields%interp_PRZ(time,i_elm,[var_rho,var_Ti,var_Te],3,st(1),st(2),phi, &
                         P,P_s,P_t,P_phi,P_time,                               &
                         R,R_s,R_t,Z,Z_s,Z_t)
  ii_Ti = 2; ii_Te = 3
  Te_jor = P(ii_Te)
  Ti_jor = P(ii_Ti)
#else
  call fields%interp_PRZ(time,i_elm,[var_rho,var_T],2,st(1),st(2),phi,         &
                         P,P_s,P_t,P_phi,P_time,                               &
                         R,R_s,R_t,Z,Z_s,Z_t)
  ii_Ti = 2; ii_Te = 2
  Te_jor = P(2)*0.5d0
  Ti_jor = P(2)*0.5d0
#endif

  ! density
  if (present(n_e_raw)) n_e_raw = P(1) * n_norm
  n_e = corr_neg_dens(P(1)) * n_norm

  ! temperatures
  ! compute Te (required)
  if (present(T_e_raw)) T_e_raw = Te_jor * T_norm
  T_e = corr_neg_temp(Te_jor*2.d0) * 0.5d0  * T_norm

  ! compute Ti only if requested (in 1-T, this is the same component anyway)
  if (need_Ti) then
    if (present(T_i_raw)) T_i_raw = Ti_jor * T_norm
    if (present(T_i))     T_i     = corr_neg_temp(Ti_jor*2.d0) * 0.5d0 * T_norm
  end if

  ! gradients (only if requested)
  if (need_grad) then
    xjac = jac(R_s,R_t,Z_s,Z_t)
    inv_xjac = 0.d0
    if (abs(xjac) > EPS) inv_xjac = 1.d0/xjac
    
    inv_R = 0.d0
    if (abs(R) > EPS) inv_R = 1.d0/R
    
    if (.not. with_TiTe) then
      ! 1-T case
      tmp_g = grad_of(ii_Te, inv_xjac, inv_R, R_s, R_t, Z_s, Z_t, P_s, P_t, P_phi, T_norm*0.5d0)
      if (present(grad_T_e)) grad_T_e = tmp_g
      if (present(grad_T_i)) grad_T_i = tmp_g
    else
      ! 2-T case
      if (present(grad_T_i)) grad_T_i = grad_of(ii_Ti, inv_xjac, inv_R, R_s, R_t, Z_s, Z_t, P_s, P_t, P_phi, T_norm)
      if (present(grad_T_e)) grad_T_e = grad_of(ii_Te, inv_xjac, inv_R, R_s, R_t, Z_s, Z_t, P_s, P_t, P_phi, T_norm)
    end if
    
    if(present(grad_raw_n_e)) grad_raw_n_e = grad_of(1, inv_xjac, inv_R, R_s, R_t, Z_s, Z_t, P_s, P_t, P_phi, n_norm)

  end if

  
  contains

    ! Helper function using suffix to avoid masking parent variables
    pure function grad_of(ii_, inv_xjac_, inv_R_, R_s_, R_t_, Z_s_, Z_t_, &
                          P_s_, P_t_, P_phi_, norm_) result(g)
      integer, intent(in)             :: ii_
      real*8, intent(in)              :: inv_xjac_, inv_R_, norm_
      real*8, intent(in)              :: R_s_, R_t_, Z_s_, Z_t_
      real*8, intent(in)              :: P_s_(:), P_t_(:), P_phi_(:)
      real*8                          :: g(3)
  
      g(1) = norm_ * ((  P_s_(ii_) * Z_t_ - P_t_(ii_) * Z_s_) * inv_xjac_)
      g(2) = norm_ * ((- P_s_(ii_) * R_t_ + P_t_(ii_) * R_s_) * inv_xjac_)
      g(3) = norm_ * (   P_phi_(ii_) * inv_R_ )
    end function grad_of

end subroutine calc_NeTeTi


!> calculate vpar in si units (note that it should still be multiplied by the norm of the B field to be si)
subroutine calc_vpar(fields, time, i_elm, st, phi, vpar)
  use phys_module, only: central_density, central_mass
  use constants
  class(fields_base), intent(in)                    :: fields
  integer, intent(in)                               :: i_elm
  real*8, intent(in)                                :: time, st(2), phi
  real*8, intent(out)                               :: vpar !< parallel velocity [m/s / T] (multiply by norm2(B) still to get [m/s])

  real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t, xjac
  real*8               :: v_norm !< vpar normalisation

  call fields%interp_PRZ(time,i_elm,[var_vpar],1,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)

  v_norm = 1.d0/sqrt(MU_ZERO*central_mass*central_density*1.d20*atomic_mass_unit)
  vpar = P(1)*v_norm

end subroutine calc_vpar

!> Calculate densities and temperature(s) for all species including ions
!> For impurities, coronal equilibrium is assumed. Note that you will need adas data to be initialized first
!> (call init_imp_adas).
!>
!> TODO make this "pure" subroutine but currently imp_cor%interp_linear prevents this.
subroutine calc_NjTj(fields, time, i_elm, st, phi, m_i_over_m_imp, ne, te, ni, ti)
  use phys_module, only: central_density, imp_cor
  use constants
  class(fields_base), intent(in) :: fields
  integer, intent(in)            :: i_elm
  real*8, intent(in)             :: time, st(2), phi
  real*8, intent(in)             :: m_i_over_m_imp !< main ion mass / mass of the impurity (used only if impurities present)
  real*8, intent(out)            :: ne             !< electron density [m^-3]
  real*8, intent(out)            :: te             !< electron temperature [K]
  real*8, intent(out)            :: ni(:)          !< ion densities [m^-3] (for each charge state: [n_main, n_imp0, n_imp+1, ...])
  real*8, intent(out)            :: ti             !< ion temperature [K]
  
  real*8, dimension(4) :: P, P_s, P_t, P_phi, P_time
  real*8               :: R, R_s, R_t, Z, Z_s, Z_t
  integer :: i

  ! Interpolate needed quantities depending on case and evaluate temperature(s)
  if(with_TiTe) then
     if(with_impurities) then
        call fields%interp_PRZ(time,i_elm,[var_rho,var_Te,var_rhoimp,var_Ti],4,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
        Ti = corr_neg_temp(P(4)*2.d0)/(2.d0*K_BOLTZ*MU_ZERO*central_density*1.d20)
     else
        call fields%interp_PRZ(time,i_elm,[var_rho,var_Te,var_Ti],3,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
        Ti = corr_neg_temp(P(3)*2.d0)/(2.d0*K_BOLTZ*MU_ZERO*central_density*1.d20)
     end if
     Te = corr_neg_temp(P(2)*2.d0)/(2.d0*K_BOLTZ*MU_ZERO*central_density*1.d20)

  else
     if(with_impurities) then
        call fields%interp_PRZ(time,i_elm,[var_rho,var_T,var_rhoimp],3,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
     else
        call fields%interp_PRZ(time,i_elm,[var_rho,var_T],2,st(1),st(2),phi,P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
     end if

     Te = corr_neg_temp(P(2))/(2.d0*K_BOLTZ*MU_ZERO*central_density*1.d20)
     Ti = Te
  end if

  if(with_impurities) then

     ni(1) = corr_neg_dens( P(1) - P(3) ) * central_density * 1d20 ! main ion density [1/m^3], capped against negative
     ne = ni(1)

     ! Assume single impurity species and store their densities to ni array as [n_main, n_imp0, n_imp+1, ...]
     call imp_cor(1)%interp_linear(density=20.d0,temperature=log10(Te),p_out=ni(2:size(ni)))

     ni(2:size(ni)) = central_density*1.d20 * m_i_over_m_imp * ni(2:size(ni)) * P(3)
     ne = ne + sum( ni(2:size(ni)) * ( (/ (i, i=0,size(ni)-2, 1) /) ) )
  else
     ne    = corr_neg_dens(P(1)) * central_density * 1d20
     ni(1) = ne
  end if

end subroutine calc_NjTj


!> Calculates the gyro-averaged electric fields from a set of particles (representing the gyro-orbit)
subroutine calc_gyro_average_E(fields, time, particles, n_phases, E_average)
  use phys_module, only: F0, mode, central_mass, central_density
  use constants, only: mu_zero, atomic_mass_unit
  use mod_particle_types
  ! Routine parameters
  class(fields_base), intent(in)             :: fields
  real*8, intent(in)                         :: time
  type(particle_kinetic_leapfrog),intent(in) :: particles(n_phases)
  integer, intent(in)                        :: n_phases     ! the number of points used in the gyro orbit average
  real*8, intent(inout)                      :: E_average(3) !< Electric field [V/m]
  ! Internal parameters
  integer, parameter :: i_var(1) = [2]
  real*8             :: P(1), P_s(1), P_t(1), P_phi(1), P_time(1) ! Placeholder for evaluating variables and derivatives locally
  real*8             :: E(3), R, R_s, R_t, Z, Z_s, Z_t
  real*8             :: inv_st_jac, R_inv
  real*8             :: U, U_R, U_Z, U_phi, t_norm
  integer            :: i

  !!!!!!!!!!!!! careful Ptime is not yet defined !!!!!!!!!!!!!!!!!!!!!!

  if (n_phases .lt. 1) return   ! return E_average as it was

  t_norm  = sqrt(mu_zero * ATOMIC_MASS_UNIT * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  ! Interpolate the fields to get psi and U at the current position (and the
  ! changes u_n - u(n-1))

  E_average = 0.d0

  do i=1, n_phases

    call fields%interp_PRZ(time, particles(i)%i_elm, i_var, 1, particles(i)%st(1), particles(i)%st(2), particles(i)%x(3), P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

    P_time(1) = 0.d0

    R_inv = 1.d0/R
    inv_st_jac = 1.d0/jac(R_s,R_t,Z_s,Z_t)

    ! Calculate the derivatives to R and Z
    U_R      = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
    U_Z      = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
    U_phi    = P_phi(1)

    U   = P(1)/t_norm

    ! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
    ! See http://jorek.eu/wiki/doku.php?id=u_phi
    E     = [-F0*U_R, -F0*U_Z, -F0*U_phi*R_inv]/t_norm
    E(3)  = E(3) - R_inv*P_time(1) ! because this is not normalized with t_norm

    E_average = E_average + E

  enddo

  E_average = E_average / real(n_phases,8)

end subroutine calc_gyro_average_E

!> Calculates the plasma velocity vector
!the sum of vpar v_ExB
! May be To do: add diamagnetic drift component.
subroutine calc_vvector(fields, time, i_elm, st, phi, vvector) 
use phys_module, only: F0, mode, central_mass, central_density
use constants, only: mu_zero, atomic_mass_unit
use mod_coordinate_transforms, only: transform_derivatives_st_to_RZ
! Routine parameters
class(fields_base), intent(in) :: fields
real*8, intent(in)  :: time
integer, intent(in) :: i_elm !< JOREK element index
real*8, intent(in)  :: st(2) !< element-local coordinates
real*8, intent(in)  :: phi !< toroidal angle
! real*8, intent(out) :: E(3) !< Electric field [V/m]
! real*8, intent(out) :: B(3) !< Magnetic field [T]
! real*8, intent(out) :: psi !< psi in JOREK units
! real*8, intent(out) :: u !< velocity stream function in m/s
real*8, intent(out) :: vvector(3) !v [v_R, v_Z, v_phi] in m/s
! Internal parameters
integer, parameter :: i_var(3) = [1,2,7]
real*8             :: P(3), P_s(3), P_t(3), P_phi(3), P_time(3) ! Placeholder for evaluating variables and derivatives locally
! Values
real*8             :: R, R_s, R_t, Z, Z_s, Z_t
! Others
real*8             :: inv_st_jac, R_inv
real*8             :: psi_R, psi_Z, U_R, U_Z, U_phi, t_norm
real*8             :: vpar, v_R, v_Z, v_phi
t_norm  = sqrt(mu_zero * ATOMIC_MASS_UNIT * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

! Interpolate the fields to get psi and U at the current position (and the
! changes u_n - u(n-1))
call fields%interp_PRZ(time, i_elm, i_var, 3, st(1), st(2), phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)

R_inv = 1.d0/R
inv_st_jac = 1.d0/jac(R_s,R_t,Z_s,Z_t)

! Calculate the derivatives to R and Z
psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
U_R      = (  P_s(2) * Z_t - P_t(2) * Z_s ) * inv_st_jac
U_Z      = (- P_s(2) * R_t + P_t(2) * R_s ) * inv_st_jac
U_phi    = P_phi(2)

! Calculate the velocity vector (see http://jorek.eu/wiki/doku.php?id=reduced_mhd)
vpar  = P(3)
v_R   = -R * U_Z + vpar * R_inv *psi_Z
v_Z   = R * U_R - vpar * R_inv *psi_R
v_phi = F0 * vpar * R_inv

vvector = [v_R, v_Z, v_phi] / t_norm

end subroutine calc_vvector


pure function rot_tmp(x,A,dA) result(rotA)
  implicit none
  real*8, intent(in)  :: x(3), A(3), dA(3,3)
  real*8              :: rotA(3)
  rotA(1) = dA(3,2) - dA(2,3) / x(1)
  rotA(2) = dA(1,3)/x(1) - dA(3,1) - A(3) / x(1)
  rotA(3) = dA(2,1) - dA(1,2)
  return
end

pure subroutine calc_RK4_analytic(fields, R, Z, phi, A_out, dA_out, B_out, dB_out, B_norm, dB_norm, bn, dBn, E)
  use phys_module, only: mode, central_mass, central_density
  use constants, only: mu_zero, atomic_mass_unit
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
  use constants, only: mu_zero, atomic_mass_unit
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

  t_norm  = sqrt(mu_zero * ATOMIC_MASS_UNIT * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

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

  RZjac    = jac(R_s,R_t,Z_s,Z_t)

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

  dBn(1) = 1.d0 /(R**2 * Bn) * (psi_R * psi_RR   + psi_Z * psi_RZ) - Bn / R
  dBn(2) = 1.d0 /(R**2 * Bn) * (psi_R * psi_RZ   + psi_Z * psi_ZZ)
  dBn(3) = 1.d0 /(R**2 * Bn) * (psi_R * psi_Rphi + psi_Z * psi_Zphi)

  dBnorm(1,1) = dB(1,1) / Bn - B(1) * dBn(1) / Bn2
  dBnorm(1,2) = dB(1,2) / Bn - B(1) * dBn(2) / Bn2
  dBnorm(1,3) = dB(1,3) / Bn - B(1) * dBn(3) / Bn2
  dBnorm(2,1) = dB(2,1) / Bn - B(2) * dBn(1) / Bn2
  dBnorm(2,2) = dB(2,2) / Bn - B(2) * dBn(2) / Bn2
  dBnorm(2,3) = dB(2,3) / Bn - B(2) * dBn(3) / Bn2
  dBnorm(3,1) = dB(3,1) / Bn - B(3) * dBn(1) / Bn2
  dBnorm(3,2) = dB(3,2) / Bn - B(3) * dBn(2) / Bn2
  dBnorm(3,3) = dB(3,3) / Bn - B(3) * dBn(3) / Bn2

! The local electric field, obtained from E=-Grad (u F0)-\partial_t A
! See http://jorek.eu/wiki/doku.php?id=u_phi
  E     = [-F0*U_R, -F0*U_Z, -F0*U_phi*R_inv]/t_norm
  E(3)  = E(3) - R_inv*P_time(1) ! because this is not normalized with t_norm

end subroutine calc_RK4

subroutine check_consistency_RK4(fields, i_elm, st)
  use phys_module, only: F0, mode, central_mass, central_density
  use constants, only: mu_zero, atomic_mass_unit
  use mod_find_rz_nearby
  class(fields_base), intent(in) :: fields
  real*8  :: time
  integer :: i_elm       !< JOREK element index
  real*8  :: st(2)       !< element-local coordinates

  real*8  :: phi         !< toroidal angle
  real*8  :: E(3)        !< Electric field [V/m]
  real*8  :: A(3)        !< vector potential [T.m]
  real*8  :: dA(3,3)     !< derivatives of vector potential [T]
  real*8  :: B(3)        !< Magnetic field [T]
  real*8  :: dB(3,3)     !< derivatives of magnetic field [T/m]
  real*8  :: Bnorm(3)    !< normalised magnetic field vector
  real*8  :: dBnorm(3,3) !< derivatives of normalised magnetic field vector [1/m]
  real*8  :: Bn          !< Magnetic field amplitude [T]
  real*8  :: dBn(3)      !< derivatives of magnetic field amplitude [T/m]

! Internal parameters
  integer, parameter :: i_var(2) = [1,2]
  real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_ss(2), P_st(2), P_tt(2), P_sphi(2), P_tphi(2)
  real*8             :: P_time(2), P_stime(2), P_ttime(2), bn2
  real*8             :: x(3), R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
! Others
  real*8             :: inv_st_jac, R_inv, RZjac, RZjac_R, RZjac_Z
  real*8             :: psi, psi_R, psi_Z, psi_RR, psi_ZZ, psi_RZ, psi_Rphi, psi_Zphi
  real*8             :: U, U_R, U_Z, U_phi, t_norm

  integer            :: i_elm_p, i_elm_m, ifail
  real*8             :: st_p(2), st_m(2), Rout, Zout, R_p, Z_p, R_m, Z_m, phi_m, phi_p, delta, error
  real*8             :: A_p(3), dA_p(3,3), B_p(3), db_p(3,3), Bnorm_p(3), dBnorm_p(3,3), Bn_p, dBn_p(3), E_p(3)
  real*8             :: A_m(3), dA_m(3,3), B_m(3), db_m(3,3), Bnorm_m(3), dBnorm_m(3,3), Bn_m, dbn_m(3), E_m(3)
  logical            :: verbose = .false.

  time = 0.d0
  phi  = 0.d0

  delta = 1.d-5

  call fields%interp_PRZ_2(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, &
                         P_time, P_ss, P_st, P_tt, P_sphi, P_tphi, P_stime, P_ttime,   &
                         R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_ss, Z_st, Z_tt)

  call calc_RK4(fields, time, i_elm, st, phi, A, dA, B, dB, Bnorm, dBnorm, bn, dBn, E)

  R_p = R + delta
  Z_p = Z
  call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_p, Z_p, st_p(1), st_p(2), i_elm_p, ifail)
  call calc_RK4(fields, time, i_elm_p, st_p, phi, A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dBn_p, E_p)

  R_m =  R - delta
  Z_m = Z
  call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_m, Z_m, st_m(1), st_m(2), i_elm_m, ifail)
  call calc_RK4(fields, time, i_elm_m, st_m, phi, A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dBn_m, E_m)

  if (verbose) then
  write(*,*) 'RK4 consistency check : R dependence'
    write(*,'(A,8e18.10)') 'A(1),  dA(1,1)  : ',A(1), dA(1,1), (A_p(1) - A_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'A(2),  dA(2,1)  : ',A(2), dA(2,1), (A_p(2) - A_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'A(3),  dA(3,1)  : ',A(3), dA(3,1), (A_p(3) - A_m(3))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'bn,    dbn(1)   : ',bn,   dbn(1),  (bn_p   - bn_m)  / (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(1),  dB(1,1)  : ',B(1), dB(1,1), (B_p(1) - B_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(2),  dB(2,1)  : ',B(2), dB(2,1), (B_p(2) - B_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(3),  dB(3,1)  : ',B(3), dB(3,1), (B_p(3) - B_m(3))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(1),  dBnorm(1,1)  : ',Bnorm(1), dBnorm(1,1), (Bnorm_p(1) - Bnorm_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(2),  dBnorm(2,1)  : ',Bnorm(2), dBnorm(2,1), (Bnorm_p(2) - Bnorm_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(3),  dBnorm(3,1)  : ',Bnorm(3), dBnorm(3,1), (Bnorm_p(3) - Bnorm_m(3))/ (2.d0*delta)
  endif

  error =         sum(abs(dA(:,1) - (A_p(:) - A_m(:))/(2.d0*delta)))
  error = error + sum(abs(dB(:,1) - (B_p(:) - B_m(:))/(2.d0*delta)))
  error = error +     abs(dbn(1)  - (bn_p   - bn_m)  /(2.d0*delta))
  error = error + sum(abs(dBnorm(:,1) - (Bnorm_p(:) - Bnorm_m(:))/(2.d0*delta)))

  R_p = R
  Z_p = Z + delta
  call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_p, Z_p, st_p(1), st_p(2), i_elm_p, ifail)
  call calc_RK4(fields, time, i_elm_p, st_p, phi, A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dBn_p, E_p)

  R_m =  R
  Z_m = Z - delta
  call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_m, Z_m, st_m(1), st_m(2), i_elm_m, ifail)
  call calc_RK4(fields, time, i_elm_m, st_m, phi, A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dBn_m, E_m)

  if (verbose) then
  write(*,*) 'RK4 consistency check : Z dependence'
    write(*,'(A,8e18.10)') 'A(1),  dA(1,2)  : ',A(1), dA(1,2), (A_p(1) - A_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'A(2),  dA(2,2)  : ',A(2), dA(2,2), (A_p(2) - A_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'A(3),  dA(3,2)  : ',A(3), dA(3,2), (A_p(3) - A_m(3))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'bn,    dbn(2)   : ',bn,   dbn(2),  (bn_p   - bn_m)  / (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(1),  dB(1,2)  : ',B(1), dB(1,2), (B_p(1) - B_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(2),  dB(2,2)  : ',B(2), dB(2,2), (B_p(2) - B_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(3),  dB(3,2)  : ',B(3), dB(3,2), (B_p(3) - B_m(3))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(1),  dBnorm(1,2)  : ',Bnorm(1), dBnorm(1,2), (Bnorm_p(1) - Bnorm_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(2),  dBnorm(2,2)  : ',Bnorm(2), dBnorm(2,2), (Bnorm_p(2) - Bnorm_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(3),  dBnorm(3,2)  : ',Bnorm(3), dBnorm(3,2), (Bnorm_p(3) - Bnorm_m(3))/ (2.d0*delta)
  endif

  error = error + sum(abs(dA(:,2) - (A_p(:) - A_m(:))/(2.d0*delta)))
  error = error + sum(abs(dB(:,2) - (B_p(:) - B_m(:))/(2.d0*delta)))
  error = error +     abs(dbn(2)  - (bn_p   - bn_m)  /(2.d0*delta))
  error = error + sum(abs(dBnorm(:,2) - (Bnorm_p(:) - Bnorm_m(:))/(2.d0*delta)))

R_p = R 
Z_p = Z
phi_p = phi + delta
call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_p, Z_p, st_p(1), st_p(2), i_elm_p, ifail)
call calc_RK4(fields, time, i_elm_p, st_p, phi_p, A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dBn_p, E_p)

R_m = R
Z_m = Z
phi_m = phi - delta
call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_m, Z_m, st_m(1), st_m(2), i_elm_m, ifail)
call calc_RK4(fields, time, i_elm_m, st_m, phi_m, A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dBn_m, E_m)

if (verbose) then
  write(*,*) 'RK4 consistency check : phi dependence'
  write(*,'(A,8e18.10)') 'A(1),  dA(1,3)  : ',A(1), dA(1,3), (A_p(1) - A_m(1))/ (2.d0*delta)
  write(*,'(A,8e18.10)') 'A(2),  dA(2,3)  : ',A(2), dA(2,3), (A_p(2) - A_m(2))/ (2.d0*delta)
  write(*,'(A,8e18.10)') 'A(3),  dA(3,3)  : ',A(3), dA(3,3), (A_p(3) - A_m(3))/ (2.d0*delta)
  write(*,'(A,8e18.10)') 'bn,    dbn(3)   : ',bn,   dbn(3),  (bn_p   - bn_m)  / (2.d0*delta)
  write(*,'(A,8e18.10)') 'B(1),  dB(1,3)  : ',B(1), dB(1,3), (B_p(1) - B_m(1))/ (2.d0*delta)
  write(*,'(A,8e18.10)') 'B(2),  dB(2,3)  : ',B(2), dB(2,3), (B_p(2) - B_m(2))/ (2.d0*delta)
  write(*,'(A,8e18.10)') 'B(3),  dB(3,3)  : ',B(3), dB(3,3), (B_p(3) - B_m(3))/ (2.d0*delta)
  write(*,'(A,8e18.10)') 'Bnorm(1),  dBnorm(1,3)  : ',Bnorm(1), dBnorm(1,3), (Bnorm_p(1) - Bnorm_m(1))/ (2.d0*delta)
  write(*,'(A,8e18.10)') 'Bnorm(2),  dBnorm(2,3)  : ',Bnorm(2), dBnorm(2,3), (Bnorm_p(2) - Bnorm_m(2))/ (2.d0*delta)
  write(*,'(A,8e18.10)') 'Bnorm(3),  dBnorm(3,3)  : ',Bnorm(3), dBnorm(3,3), (Bnorm_p(3) - Bnorm_m(3))/ (2.d0*delta)
endif

error =         sum(abs(dA(:,3) - (A_p(:) - A_m(:))/(2.d0*delta)))
error = error + sum(abs(dB(:,3) - (B_p(:) - B_m(:))/(2.d0*delta)))
error = error +     abs(dbn(3)  - (bn_p   - bn_m)  /(2.d0*delta))
error = error + sum(abs(dBnorm(:,3) - (Bnorm_p(:) - Bnorm_m(:))/(2.d0*delta)))

  write(*,*) 'RK4 consistency : error : ',error

  return
end subroutine check_consistency_RK4

pure subroutine calc_Qin_analytic(fields, R, Z, phi, A_out, dA_out, B_out, dB_out, B_norm, dB_norm, bn, dBn, E)
  use phys_module, only: mode, central_mass, central_density
  use constants, only: mu_zero, atomic_mass_unit
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
  use constants, only: mu_zero, atomic_mass_unit
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

  t_norm  = sqrt(mu_zero * ATOMIC_MASS_UNIT * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  call fields%interp_PRZ_2(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, &
                         P_time, P_ss, P_st, P_tt, P_sphi, P_tphi, P_stime, P_ttime,   &
                         R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_ss, Z_st, Z_tt)

  R_inv = 1.d0/R
  inv_st_jac = 1.d0/jac(R_s,R_t,Z_s,Z_t)

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

  RZjac    = jac(R_s,R_t,Z_s,Z_t)

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

  dBn(1) = 1.d0 /(R**2 * Bn) * (psi_R * psi_RR   + psi_Z * psi_RZ) - Bn / R
  dBn(2) = 1.d0 /(R**2 * Bn) * (psi_R * psi_RZ   + psi_Z * psi_ZZ)
  dBn(3) = 1.d0 /(R**2 * Bn) * (psi_R * psi_Rphi + psi_Z * psi_Zphi)

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

subroutine check_consistency_Qin(fields, i_elm, st)
  use phys_module, only: F0, mode, central_mass, central_density
  use constants, only: mu_zero, atomic_mass_unit
  use mod_find_rz_nearby
  class(fields_base), intent(in) :: fields
  real*8  :: time
  integer :: i_elm       !< JOREK element index
  real*8  :: st(2)       !< element-local coordinates

  real*8  :: phi         !< toroidal angle
  real*8  :: E(3)        !< Electric field [V/m]
  real*8  :: A(3)        !< vector potential [T.m]
  real*8  :: dA(3,3)     !< derivatives of vector potential [T]
  real*8  :: B(3)        !< Magnetic field [T]
  real*8  :: dB(3,3)     !< derivatives of magnetic field [T/m]
  real*8  :: Bnorm(3)    !< normalised magnetic field vector
  real*8  :: dBnorm(3,3) !< derivatives of normalised magnetic field vector [1/m]
  real*8  :: Bn          !< Magnetic field amplitude [T]
  real*8  :: dBn(3)      !< derivatives of magnetic field amplitude [T/m]

  ! Internal parameters
  integer, parameter :: i_var(2) = [1,2]
  real*8             :: P(2), P_s(2), P_t(2), P_phi(2), P_ss(2), P_st(2), P_tt(2), P_sphi(2), P_tphi(2)
  real*8             :: P_time(2), P_stime(2), P_ttime(2), bn2
  real*8             :: x(3), R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
  ! Others
  real*8             :: inv_st_jac, R_inv, RZjac, RZjac_R, RZjac_Z
  real*8             :: psi, psi_R, psi_Z, psi_RR, psi_ZZ, psi_RZ, psi_Rphi, psi_Zphi
  real*8             :: U, U_R, U_Z, U_phi, t_norm

  integer            :: i_elm_p, i_elm_m, ifail
  real*8             :: st_p(2), st_m(2), Rout, Zout, R_p, Z_p, R_m, Z_m, delta, error
  real*8             :: A_p(3), dA_p(3,3), B_p(3), db_p(3,3), Bnorm_p(3), dBnorm_p(3,3), Bn_p, dBn_p(3), E_p(3)
  real*8             :: A_m(3), dA_m(3,3), B_m(3), db_m(3,3), Bnorm_m(3), dBnorm_m(3,3), Bn_m, dbn_m(3), E_m(3)
  logical            :: verbose = .false.

  time = 0.d0
  phi  = 0.d0

  delta = 1.d-5

  call fields%interp_PRZ_2(time, i_elm, i_var, 2, st(1), st(2), phi, P, P_s, P_t, P_phi, &
                          P_time, P_ss, P_st, P_tt, P_sphi, P_tphi, P_stime, P_ttime,   &
                          R, R_s, R_t, R_ss, R_st, R_tt, Z, Z_s, Z_t, Z_ss, Z_st, Z_tt)

  call calc_Qin(fields, time, i_elm, st, phi, A, dA, B, dB, Bnorm, dBnorm, bn, dBn, E)

  R_p = R + delta
  Z_p = Z
  call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_p, Z_p, st_p(1), st_p(2), i_elm_p, ifail)
  call calc_Qin(fields, time, i_elm_p, st_p, phi, A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dBn_p, E_p)

  R_m =  R - delta
  Z_m = Z
  call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_m, Z_m, st_m(1), st_m(2), i_elm_m, ifail)
  call calc_Qin(fields, time, i_elm_m, st_m, phi, A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dBn_m, E_m)

  if (verbose) then
    write(*,*) 'Qin consistency check : '
    write(*,'(A,8e18.10)') 'A(1),  dA(1,1)  : ',A(1), dA(1,1), (A_p(1) - A_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'A(2),  dA(2,1)  : ',A(2), dA(2,1), (A_p(2) - A_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'A(3),  dA(3,1)  : ',A(3), dA(3,1), (A_p(3) - A_m(3))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'bn,    dbn(1)   : ',bn,   dbn(1),  (bn_p   - bn_m)  / (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(1),  dB(1,1)  : ',B(1), dB(1,1), (B_p(1) - B_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(2),  dB(2,1)  : ',B(2), dB(2,1), (B_p(2) - B_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(3),  dB(3,1)  : ',B(3), dB(3,1), (B_p(3) - B_m(3))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(1),  dBnorm(1,1)  : ',Bnorm(1), dBnorm(1,1), (Bnorm_p(1) - Bnorm_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(2),  dBnorm(2,1)  : ',Bnorm(2), dBnorm(2,1), (Bnorm_p(2) - Bnorm_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(3),  dBnorm(3,1)  : ',Bnorm(3), dBnorm(3,1), (Bnorm_p(3) - Bnorm_m(3))/ (2.d0*delta)
  endif

  error =         sum(abs(dA(:,1) - (A_p(:) - A_m(:))/(2.d0*delta)))
  error = error + sum(abs(dB(:,1) - (B_p(:) - B_m(:))/(2.d0*delta)))
  error = error +     abs(dbn(1)  - (bn_p   - bn_m)  /(2.d0*delta))
  error = error + sum(abs(dBnorm(:,1) - (Bnorm_p(:) - Bnorm_m(:))/(2.d0*delta)))

  R_p = R
  Z_p = Z + delta
  call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_p, Z_p, st_p(1), st_p(2), i_elm_p, ifail)
  call calc_Qin(fields, time, i_elm_p, st_p, phi, A_p, dA_p, B_p, dB_p, Bnorm_p, dBnorm_p, bn_p, dBn_p, E_p)

  R_m =  R
  Z_m = Z - delta
  call find_RZ_nearby(fields%node_list, fields%element_list, R, Z, st(1), st(2), i_elm, R_m, Z_m, st_m(1), st_m(2), i_elm_m, ifail)
  call calc_Qin(fields, time, i_elm_m, st_m, phi, A_m, dA_m, B_m, dB_m, Bnorm_m, dBnorm_m, bn_m, dBn_m, E_m)

  if (verbose) then
    write(*,'(A,8e18.10)') 'A(1),  dA(1,2)  : ',A(1), dA(1,2), (A_p(1) - A_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'A(2),  dA(2,2)  : ',A(2), dA(2,2), (A_p(2) - A_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'A(3),  dA(3,2)  : ',A(3), dA(3,2), (A_p(3) - A_m(3))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'bn,    dbn(2)   : ',bn,   dbn(2),  (bn_p   - bn_m)  / (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(1),  dB(1,2)  : ',B(1), dB(1,2), (B_p(1) - B_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(2),  dB(2,2)  : ',B(2), dB(2,2), (B_p(2) - B_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'B(3),  dB(3,2)  : ',B(3), dB(3,2), (B_p(3) - B_m(3))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(1),  dBnorm(1,2)  : ',Bnorm(1), dBnorm(1,2), (Bnorm_p(1) - Bnorm_m(1))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(2),  dBnorm(2,2)  : ',Bnorm(2), dBnorm(2,2), (Bnorm_p(2) - Bnorm_m(2))/ (2.d0*delta)
    write(*,'(A,8e18.10)') 'Bnorm(3),  dBnorm(3,2)  : ',Bnorm(3), dBnorm(3,2), (Bnorm_p(3) - Bnorm_m(3))/ (2.d0*delta)
  endif

  error = error + sum(abs(dA(:,2) - (A_p(:) - A_m(:))/(2.d0*delta)))
  error = error + sum(abs(dB(:,2) - (B_p(:) - B_m(:))/(2.d0*delta)))
  error = error +     abs(dbn(2)  - (bn_p   - bn_m)  /(2.d0*delta))
  error = error + sum(abs(dBnorm(:,2) - (Bnorm_p(:) - Bnorm_m(:))/(2.d0*delta)))

  write(*,*) 'Qin consistency : error : ',error

  return
end subroutine check_consistency_Qin

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
  use constants, only: mu_zero,atomic_mass_unit
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
  U = F0*U/sqrt(mu_zero*ATOMIC_MASS_UNIT*central_mass*central_density*1.d20)

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

!> calculates the jacobian R_s*Z_t - R_t*Z_s
!> returns small number if jac = 0 (to avoid NaNs, but it is of course not correct)
real*8 function jac(R_s,R_t,Z_s,Z_t)
  implicit none
  real*8, intent(in)           :: R_s,R_t,Z_s,Z_t

  real*8, parameter :: tol = 1.d-20

  jac = (R_s * Z_t - R_t * Z_s)
  
  !> WARNING: this is to avoid a NaN from ruining the simulation, but it is not correct!
  if (abs(jac) < tol) then
    !$omp critical
    write(*,"(A)") "ERROR: jacobian=0 in particles/mod_fields.f90. Returning small number instead to avoid NaNs"
    write(*,"(A)") "This should normally not happen. Likely something was sampled exactly at the grid_axis."
    write(*,"(A)") "Please find and solve the underlying problem."
    !$omp end critical
    jac = sign(tol,jac)
  endif
end function jac

end module mod_fields
