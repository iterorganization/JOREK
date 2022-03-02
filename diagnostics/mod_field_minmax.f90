!> mod_field_minmax contains variables and procedures for 
!> for finding the local and global maximum and minimum
!> values of node list fields
module mod_fields_minmax
use mod_rootfinding, only: fun_dfun_3d
implicit none
private
public ::

!> Variables and datatypes -------------------------------
!> wrap the interpolation function in a class for passing
!> it as an argument to the newton root finding method
type,extends(fun_dfun_3d_loc) :: fun_interp_PRZ
  use data_structure, only: type_node_list
  use data_structure, only: type_element_list
  implicit none
  type(type_node_list)    :: node_list
  type(type_element_list) :: element_list
  contains
  procedure,pass :: f_df => f_df_interp_PRZ
end type fun_interp_PRZ

!> Interfaces---------------------------------------------
!> interface foro fun_interp_PRZ constructor
interface fun_interp_PRZ
  module procedure init_fun_interp_PRZ
end interface fun_interp_PRZ
contains
!> Procedures --------------------------------------------
!> generate random mesh in the s,t,phi coordinates
!> inputs:
!>   n_trials: (integer) number of mesh elements
!>   rngs:     (type_rng)(n_threads) random number generators
!> outputs:
!>   mesh:     (real8)(3,n_trials) s,t,phi mesh
subroutine generate_random_stphi_mesh(n_trials,rngs,mesh)
  use constants, only: TWOPI
  use mod_rng,   only: type_rng
  !$ use omp_lib
  implicit none
  !> inputs:
  integer,intent(in)                                  :: n_trials
  class(type_rng),dimension(:),allocatable,intent(in) :: rngs
  !> outputs:
  real*8,dimension(3,n_trials),intent(out)            :: mesh
  !> variables
  integer :: ii,thread_id
  !> generate mesh
  !$omp parallel default(private) firstprivate(n_trials) &
  !$omp shared(rngs,mesh)
  thread_id = 1;
  !$ thread_id = omp_get_thread_num()
  !$omp do
  do ii=1,n_trials
    call rngs(thread_id)%next(mesh(:,ii))
    mesh(3,ii) = TWOPI*mesh(3,ii)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine generate_random_stphi_mesh

!> generate equidistant mesh in the s,t,phi coordinates
!> inputs:
!>   n_s:   (integer) number of points in the s local coords.
!>   n_t:   (integer) number of points in the t local coords.
!>   n_phi: (integer) number of toroidal angle points
!> outputs:
!>   mesh:  (real8)(3,n_s*n_t*n_phi) equidistant mesh in s,t,phi
subroutine generate_equidistant_stphi_mesh(n_s,n_t,n_phi,mesh)
  use constants, only: TWOPI
  implicit none
  !> inputs:
  integer,intent(in) :: n_s,n_t,n_phi
  !> outputs:
  real*8,dimension(3,n_s*n_t*n_phi),intent(out) :: mesh
  !> variables:
  integer :: ii,jj,kk
  real*8 :: d_s,d_t,d_phi
  !> initialisation
  d_s   = 1.d0/real(n_s-1,kind=8);
  d_t   = 1.d0/real(n_t-1,kind=8);
  d_phi = TWOPI/real(n_phi-1,kind=8);
  !> generate equidistant points
  !$omp parallel do default(private) firstprivate(n_s,n_t,n_phi,&
  !$omp d_s,d_t,d_phi) shared(mesh) collapse(2)
  do ii=1,n_phi
    do jj=1,n_t
      !$omp simd
      do kk=1,n_s
        mesh(:,((ii-1)*n_t+jj-1)*n_s+kk) = (/real(kk-1,kind=8)*d_s,&
        real(jj-1,kind=8)*d_t,real(ii-1,kind=8)*d_phi/)
      enddo
      !$omp end simd
    enddo
  enddo
  !$omp end parallel do

end subroutine generate_equidistant_stphi_mesh

!> find local and global minimum and maximum values of a jorek field
!> inputs:
!> outputs:
subroutine field_minmax(field_id,n_mesh,mesh,f_interp_PRZ,minmax_list,minmax_global)
  implicit none
  !> inputs:
  class(fun_interp_PRZ),intent(in)      :: f_interp_PRZ
  integer,intent(in)                    :: fieald_id,n_mesh
  real*8,dimension(3,n_mesh),intent(in) :: mesh
  !> outputs:
  real*8,dimension(2),intent(out)       :: minmax_global
  real*8,dimension(f_interp_PRZ%element_list%n_elements),intent(out) :: min_list
  real*8,dimension(f_interp_PRZ%element_list%n_elements),intent(out) :: max_list
  !> variables:
  integer               :: ii,jj,ierr,maxit,n_elements
  real*8                :: tol
  real*8,dimension(3)   :: x_extrema,values
  real*8,dimension(3,3) :: Jac
  !> initialisation
  tol=5.d-15; maxit=10000; n_elements=f_interp_PRZ%element_list%n_elements
  minmax_list = (/1.d21,-1.d21/); min_list = 1.d21; max_list = -1.d21;
  !> find minimum and maximum
  !$omp parallel do default(private) firstprivate(n_elements,n_mesh) &
  !$omp shared(f_interp_PRZ,min_list,max_list) reduction(min:min_list) &
  !$omp reduction(max:max_list) collapse(2)
  do ii=1,n_elements
    do jj=1,n_mesh
      !> find extrema
      ierr = 0;
      call newtons_method(f_interp_PRZ,(/0.d0,0.d0,0.d0/),mesh(:,jj),&
      x_extrema,2,(/field_id,ii/),ierr,tol,maxit)
      if(ierr.ne.0) cycle
      !> compute the jacobian
      call f_interp_PRZ%f_df(values,Jac,x_extrema,2,(/field_id,ii/))
      !> check if minimum or maximum 
      if()
      if()
    enddo
  enddo
  !$omp end parallel do
end subroutine field_minmax

!> initialise the fun_interp_PRZ class
!> inputs:
!>   node_list:    (type_node_list) jorek node list
!>   element_list: (type_element_list) jorek element list
!> outputs:
!>   this: (fun_interp_PRZ) initialised fun_interp_PRZ class
function init_fun_interp_PRZ(node_list,element_list) result(this)
  use data_structure, only: type_node_list
  use data_structure, only: type_element_list
  implicit none
  !> inputs:
  type(type_node_list),intent(in)    :: node_list
  type(type_element_list).intent(in) :: element_list
  !> output:
  class(fun_interp_PRZ) :: this
  !> initialise class
  this%node_list = node_list; this%element_list = element_list;
end function init_fun_interp_PRZ

!> function used for computing the derivatives and the jacobian
!> of the derivatives of a jorek field
!> inputs:
!>   this:         (fun_interp_PRZ) fun_interp_PRZ class
!>   x:            (real8)(3) (s,t,phi) coordinates
!>   n_int_coords: (integer) must be 2
!>   int_coords:   (integer)(n_int_coords) integer coordinates:
!>                 1: element number
!>                 2: field number
!> outputs:
!>   this: (fun_interp_PRZ) fun_interp_PRZ class
!>   f:    (real8)(3) jorek field derivatives
!>   J:    (real8)(3,3) jorek field hessian
subroutine f_df_interp_PRZ(this,f,J,x,n_int_coords,int_coords)
  use mod_interp, only: interp_PRZ
  implicit none
  !> inputs-outputs:
  class(fun_interp_PRZ),intent(inout) :: this
  !> inputs:
  integer,intent(in)                         :: n_int_coords
  integer,dimension(n_int_coords),intent(in) :: int_coords
  real*8,dimension(3),intent(in)             :: x
  !> outputs:
  real*8,dimension(3),intent(out)   :: f
  real*8,dimension(3,3),intent(out) :: J
  !> variables:
  real*8 :: R,R_s,R_t,R_st,R_ss,R_tt
  real*8 :: Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
  real*8,dimension(1) :: P,P_s,P_t,P_phi
  real*8,dimension(1) :: P_st,P_ss,P_tt,P_sphi,P_tphi,P_phiphi

  !> interpolate the JOREK fields
  call interp_PRZ(this%node_list,this%element_list,int_coords(1),&
  int_coords(2),1,x(1),x(2),x(3),P,P_s,P_t,P_phi,P_st,P_ss,P_tt,&
  P_sphi,P_tphi,P_phiphi,R,R_s,R_t,R_st,R_ss,R_tt,&
  Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
  !> fill up the values and the jacobian with the first
  !> and second order derivatives
  f = (/P_s,P_t,P_phi/)
  J(:,1) = (/P_ss,P_st,P_sphi/)
  J(:,2) = (/P_st,P_tt,P_tphi/)
  J(:,3) = (/P_sphi,P_tphi,P_phiphi/)
end subroutine f_df_interp_PRZ

!>--------------------------------------------------------
end module mod_fileds_minmax
