!> mod_field_minmax contains variables and procedures for 
!> for finding the local and global maximum and minimum
!> values of node list fields
module mod_fields_minmax
use data_structure,  only: type_node_list
use data_structure,  only: type_element_list
use mod_rootfinding, only: fun_dfun_3d
implicit none
private
public :: fun_interp_PRZ
public :: generate_stphi_mesh
public :: field_minmax

!> Variables and datatypes -------------------------------
!> wrap the interpolation function in a class for passing
!> it as an argument to the newton root finding method
type,extends(fun_dfun_3d) :: fun_interp_PRZ
  type(type_node_list)    :: node_list
  type(type_element_list) :: element_list
  real*8                  :: val
  contains
  procedure,pass :: f_df => f_df_interp_PRZ
end type fun_interp_PRZ

!> Interfaces---------------------------------------------
!> interface for fun_interp_PRZ constructor
interface fun_interp_PRZ
  module procedure init_fun_interp_PRZ
end interface fun_interp_PRZ

!> interface for the mesh generation method
interface generate_stphi_mesh
  module procedure generate_random_stphi_mesh
  module procedure generate_equidistant_stphi_mesh
end interface generate_stphi_mesh
contains
!> Procedures --------------------------------------------
!> generate random mesh in the s,t,phi coordinates
!> inputs:
!>   n_trials: (integer) number of mesh elements
!>   rngs:     (type_rng)(n_threads) random number generators
!> outputs:
!>   rngs:     (type_rng)(n_threads) random number generators
!>   mesh:     (real8)(3,n_trials) s,t,phi mesh
subroutine generate_random_stphi_mesh(n_trials,rngs,mesh)
  use constants, only: TWOPI
  use mod_rng,   only: type_rng
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  class(type_rng),dimension(:),allocatable,intent(inout) :: rngs
  !> inputs:
  integer,intent(in)                                     :: n_trials
  !> outputs:
  real*8,dimension(3,n_trials),intent(out)               :: mesh
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
!>   field_id:      (integer) index of the jorek field
!>   n_mesh:        (integer) number of mesh element
!>   mesh:          (real8)(3,n_mesh) mesh element list (s,t,phi)
!>   f_interp_PRZ:  (fun_interp_PRZ) interpolation function for root finding
!> outputs:
!>   f_interp_PRZ:  (fun_interp_PRZ) interpolation function for root finding
!>   min_list:      (real8)(n_elements) minimum of each element
!>   max_list:      (real8)(n_elements) maximum of each element
!>   minmax_global: (real8)(2) global minimum and maximum of the field
subroutine field_minmax(field_id,n_mesh,mesh,f_interp_PRZ,min_list,max_list,minmax_global)
  use mod_rootfinding,    only: newtons_method
  use mod_math_operators, only: compute_eigenvalues
  implicit none
  !> inputs-outputs:
  type(fun_interp_PRZ),intent(inout)    :: f_interp_PRZ
  !> inputs:
  integer,intent(in)                    :: field_id,n_mesh
  real*8,dimension(3,n_mesh),intent(in) :: mesh
  !> outputs:
  real*8,dimension(2),intent(out)       :: minmax_global
  real*8,dimension(f_interp_PRZ%element_list%n_elements),intent(out) :: min_list
  real*8,dimension(f_interp_PRZ%element_list%n_elements),intent(out) :: max_list
  !> variables:
  integer                 :: ii,jj,ierr,maxit,n_elements
  real*8                  :: tol
  real*8,dimension(3)     :: x_extrema,values
  real*8,dimension(3,3)   :: Jac
  complex*16,dimension(3) :: eigv
  !> initialisation
  tol=5.d-16; maxit=10000; n_elements=f_interp_PRZ%element_list%n_elements
  minmax_global = (/1.d21,-1.d21/); min_list = 1.d21; max_list = -1.d21;
  !> find minimum and maximum
 ! !$omp parallel do default(private) firstprivate(n_elements,n_mesh,maxit,tol,field_id) &
 ! !$omp shared(f_interp_PRZ,mesh) reduction(min:min_list) &
 ! !$omp reduction(max:max_list) collapse(2)
  do ii=1,n_elements
    do jj=1,n_mesh
      !> find extrema
      call newtons_method(f_interp_PRZ,(/0.d0,0.d0,0.d0/),mesh(:,jj),&
      x_extrema,2,(/ii,field_id/),ierr,tol,maxit)
      call f_interp_PRZ%f_df(values,Jac,x_extrema,2,(/ii,field_id/))
      if(ierr.ne.0) cycle
      !> compute the jacobian
      call f_interp_PRZ%f_df(values,Jac,x_extrema,2,(/ii,field_id/))
      !> compute the eigenvalues of the jacobian
      call compute_eigenvalues(Jac,eigv)
      !> check for local maxima and minima: the eigenvalues must be real due
      !> to the symmetry of the hessian matrix
      if(all(real(eigv).gt.0.d0).and.all(aimag(eigv).eq.0.d0)) then
        min_list(ii) = min(min_list(ii),f_interp_PRZ%val) !< found a minimum
      elseif(all(real(eigv).lt.0.d0).and.all(aimag(eigv).eq.0.d0)) then
        max_list(ii) = max(max_list(ii),f_interp_PRZ%val) !< found a maximum
      endif
    enddo
  enddo
 ! !$omp end parallel do
  !> extract the approximate global minimum and maximum
  minmax_global = (/minval(min_list),maxval(max_list)/)
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
  type(type_element_list),intent(in) :: element_list
  !> output:
  type(fun_interp_PRZ) :: this
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
pure subroutine f_df_interp_PRZ(this,f,J,x,n_int_coords,int_coords)
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
  (/int_coords(2)/),1,x(1),x(2),x(3),P,P_s,P_t,P_phi,P_st,P_ss,P_tt,&
  P_sphi,P_tphi,P_phiphi,R,R_s,R_t,R_st,R_ss,R_tt,&
  Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
  !> fill up the values and the jacobian with the first
  !> and second order derivatives and store the value
  this%val = P(1); f = (/P_s,P_t,P_phi/)
  J(:,1) = (/P_ss,P_st,P_sphi/)
  J(:,2) = (/P_st,P_tt,P_tphi/)
  J(:,3) = (/P_sphi,P_tphi,P_phiphi/)
end subroutine f_df_interp_PRZ

!>--------------------------------------------------------
end module mod_fields_minmax
