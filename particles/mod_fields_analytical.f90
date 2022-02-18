!> fields_analytical extends fields_base implementing
!> analytical methods for computing E and B fields.
!> dummy functions for interp_PRZ and interp_PRZ_2
module mod_fields_analytical
use mod_fields, only: fields_base
implicit none
private
public :: fields_analytical

!> Variables and datatypes ----------------------------
type,extends(fields_base) :: fields_analytical
  !> plasma parameters at the magnetic axis
  real*8              :: B0,U0
  real*8,dimension(2) :: RZ0
  contains
  procedure,public :: interp_PRZ   => interp_PRZ_dummy
  procedure,public :: interp_PRZ_2 => interp_PRZ_2_dummy
  procedure,public :: calc_EBPsiU  => calc_EBPsiU_parabolic
  procedure,public :: calc_EBNormBGradBCurlbDbdt => calc_EBNormBGradBCurlbDbdt_parabolic
  procedure,public :: init_fields  => init_fields_analytical 
  procedure,public :: deallocate_fields => deallocate_fields_analytical
end type fields_analytical
!> Interfaces -----------------------------------------
contains
!> Procedures -----------------------------------------
!> Initialise the fields analytical object
!> inputs:
!>   fields:       (fields_analytical) fields to be initialised
!>   n_int_param:  (integer) N# integer parameters, must be 0
!>   n_real_param: (integer) N# real parameters, must be 4
!>   int_param:    (integer)(n_int_param) integer parameters
!>   real_param:   (real8)(n_real_param) real parameters:
!>                 1-2: R,Z position of the magnetic axis
!>                 3:   toroidal magnetic field at the magnetic axis
!>                 4:   reference electric potential intensity
!> outputs:
!>   fields:       (fields_analytical) initialised fields
subroutine init_fields_analytical(fields,n_int_param,&
n_real_param,int_param,real_param)
  implicit none
  !> inputs-outputs:
  class(fields_analytical),intent(inout) :: fields
  !> inputs:
  integer :: n_int_param,n_real_param
  integer,dimension(n_int_param),intent(in) :: int_param
  real*8,dimension(n_real_param),intent(in) :: real_param
  !> check if the number of parameters are correct and initialise
  fields%RZ0 = (/3.d0,0.d0/); fields%B0 = 3.d0; fields%U0 = 0.d0;
  if((n_int_param.eq.0).and.(n_real_param.eq.4)) then
    write(*,'(/A)') "Error initialise fields analytical: size input parameters mismatch!"
    write(*,'(/A)') "use default parameters"
    return
  endif
  fields%RZ0 = real_param(1:2); fields%B0 = real_param(3); 
  fields%U0 = real_param(4)
end subroutine init_fields_analytical

!> clean-up fields analytical
subroutine deallocate_fields_analytical(fields)
  implicit none
  class(fields_analytical),intent(inout) :: fields
  !> cleanup all parameperts
  fields%RZ0=0.d0; fields%B0=0.d0; fields%U0=0.d0;
end subroutine deallocate_fields_analytical

!> Subroutine to ocompute analytical magnetic and electric fields
!> for testing integrators. !> electric field is toroidal 
!> with intensity U*R0/R while a tokamak-like magnetic field with 
!> a poloidal flux of
!> 0.5*B0*((R-R0)**2+(Z-Z0)**2) is used.
!> inputs:
!>   fields: (fields_analytical) particle fields object
!>   time:   (real8) fields time: NOT USED
!>   i_elm:  (integer) mesh element index: NOT USED
!>   st:     (real8)(2) particle major radius and vertical position
!>   phi:    (real8) particle toroidal angle: NOT USED
!> outputs:
!>   B:   (real8)(3) magnetic field
!>   E:   (real8)(3) electric field
!>   psi: (real8) poloidal flux
!>   U:   (real8) electric potential intensity
pure subroutine calc_EBPsiU_parabolic(fields,time,i_elm,&
st,phi,E,B,psi,U)
  implicit none
  !> delcare input variables
  class(fields_analytical), intent(in) :: fields
  integer,intent(in)                   :: i_elm
  real*8,intent(in)                    :: time,phi
  real*8, dimension(2), intent(in)     :: st
  !> declare output variables:
  real*8, intent(out)               :: psi,U
  real*8, dimension(3), intent(out) :: E,B
  !> computing magnetic field
  B = fields%B0*[st(2)-fields%RZ0(2),fields%RZ0(1)-st(1),fields%RZ0(1)]/st(1)
  !> computing electric field
  E = fields%U0*[0.d0,0.d0,fields%RZ0(1)/st(1)]
  !> compute psi
  psi = 0.5*fields%B0*(dot_product(st-fields%RZ0,st-fields%RZ0))
  !> compute U
  U = fields%U0
end subroutine calc_EBPsiU_parabolic

!> The procedure computes analytical guiding ceneter
!> fields for a static electromagnetic field. The
!> electric field is toroidal with intensity U*R0/R
!> while a tokamak-like magnetic field with a poloidal flux of:
!> psi = 0.5*B0*((R-R0)**2 + (Z-Z0)**2) is used.
!> inputs:
!>   fields: (fields_analytical) particle fields object
!>   time:   (real8) fields time: NOT USED
!>   i_elm:  (integer) mesh element index: NOT USED
!>   st:     (real8)(2) particle major radius and vertical position
!>   phi:    (real8) particle toroidal angle: NOT USED
!> outputs:
!>   E:     (real8)(3) electric field
!>   b:     (real8)(3) magnetic field direction
!>   normB: (real8) magnetic intensity
!>   gradB: (real8)(3) gradient of the magnetic intensity
!>   curlb: (real8)(3) curl of the magnetic direction
!>   dbdt:  (real8)(3) magnetic direction time variation
pure subroutine calc_EBNormBGradBCurlbDbdt_parabolic(fields,time, &
  i_elm,st,phi,E,b,normB,gradB,curlb,dbdt)
  use mod_math_operators, only: cross_product
  implicit none
  !> input variables
  class(fields_analytical),intent(in) :: fields
  integer,intent(in)                  :: i_elm
  real*8,intent(in)                   :: time,phi
  real*8,dimension(2),intent(in)      :: st
  !> output variables
  real*8,intent(out)              :: normB
  real*8,dimension(3),intent(out) :: E,b,gradB,curlb,dbdt
  !> compute electric field
  E = fields%U0*[0.d0,0.d0,fields%RZ0(1)/st(1)]
  !> compute magnetic field
  b = fields%B0*[st(2)-fields%RZ0(2),fields%RZ0(1)-st(1),fields%RZ0(1)]/st(1)
  !> compute norm of the magnetic field
  normB = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
  !> compute gradient of the magnetic field
  gradB = [fields%B0*fields%B0*(st(1)-fields%RZ0(1))-normB*normB*st(1), &
          fields%B0*fields%B0*(st(2)-fields%RZ0(2)),0.d0]/(normB*st(1)*st(1))
  !> compute the magetic direction
  b = b/normB
  !> compute the curl of the magnetic field directon
  curlb = (cross_product(b,gradB) -                 &
    [0.d0,0.d0,(st(1)+fields%RZ0(1))/(st(1)*st(1))])/normB
  !> compute magnetic field time derivative
  dbdt = [0.d0,0.d0,0.d0]
end subroutine calc_EBNormBGradBCurlbDbdt_parabolic

!> dummy interpolation function returning zero
pure subroutine interp_PRZ_dummy(this,time,i_elm,i_v,n_v,s,t,phi,&
P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
  implicit none
  !> inputs:
  class(fields_analytical),intent(in) :: this
  integer,intent(in)                  :: i_elm,n_v
  integer,dimension(n_v),intent(in)   :: i_v
  real*8,intent(in)                   :: time,s,t,phi
  !> outputs:
  real*8,intent(out)                  :: R,R_s,R_t,Z,Z_s,Z_t
  real*8,dimension(n_v),intent(out)   :: P,P_s,P_t,P_time,P_phi
  !> return 0:
  R=0.d0; R_s=0.d0; R_t=0.d0; Z=0.d0; Z_s=0.d0; Z_t=0.d0;
  P=0.d0; P_s=0.d0; P_t=0.d0; P_time=0.d0; P_phi=0.d0;
end subroutine interp_PRZ_dummy

!> dummy interpolation function with second derivatives returning zero
pure subroutine interp_PRZ_2_dummy(this,time,i_elm,i_v,n_v,s,t,phi,&
P,P_s,P_t,P_phi,P_time,P_ss,P_st,P_tt,P_sphi,P_tphi,P_stime,P_ttime,&
R,R_s,R_t,R_ss,R_st,R_tt,Z,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  implicit none
  !> inputs:
  class(fields_analytical),intent(in) :: this
  integer,intent(in)                  :: i_elm,n_v
  integer,dimension(n_v),intent(in)   :: i_v
  real*8,intent(in)                   :: time,s,t,phi
  !> outputs:
  real*8,intent(out)                  :: R,R_s,R_t,R_ss,R_st,R_tt
  real*8,intent(out)                  :: Z,Z_s,Z_t,Z_ss,Z_st,Z_tt
  real*8,dimension(n_v),intent(out)   :: P,P_s,P_t,P_phi,P_time
  real*8,dimension(n_v),intent(out)   :: P_ss,P_st,P_tt,P_sphi,P_tphi
  real*8,dimension(n_v),intent(out)   :: P_stime,P_ttime
  !> return 0
  R=0.d0;R_s=0.d0;R_t=0.d0;R_ss=0.d0;R_st=0.d0;R_tt=0.d0;
  Z=0.d0;Z_s=0.d0;Z_t=0.d0;Z_ss=0.d0;Z_st=0.d0;Z_tt=0.d0;
  P=0.d0;P_s=0.d0;P_t=0.d0;P_phi=0.d0;P_time=0.d0;
  P_ss=0.d0;P_st=0.d0;P_tt=0.d0;P_sphi=0.d0;P_tphi=0.d0;
  P_stime=0.d0;P_ttime=0.d0;
end subroutine interp_PRZ_2_dummy
!>-----------------------------------------------------
end module mod_fields_analytical

