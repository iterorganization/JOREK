!> mod_field_minmax contains variables and procedures for 
!> for finding the local and global maximum and minimum
!> values of node list fields
module mod_fields_minmax
use data_structure,  only: type_node_list
use data_structure,  only: type_element_list
use mod_rootfinding, only: fun_dfun_2d,fun_dfun_3d
implicit none
private
public :: fun_interp_PRZ
public :: generate_mesh_1d
public :: field_minmax

!> Variables and datatypes -------------------------------
!> wrap the interpolation function in a class for passing
!> it as an argument to the newton root finding method
!> for axisymmetric fields only!
type,extends(fun_dfun_2d) :: fun_interp_PRZ_axisym
  type(type_node_list)    :: node_list
  type(type_element_list) :: element_list
  real*8                  :: val
  contains
  procedure,pass :: f_df => f_df_interp_PRZ_axisym
end type fun_interp_PRZ_axisym

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
!> interface for fun_interp_PRZ_axisym constructor
interface fun_interp_PRZ_axisym
  module procedure init_fun_interp_PRZ_axisym
end interface fun_interp_PRZ_axisym
!> interface for fun_interp_PRZ constructor
interface fun_interp_PRZ
  module procedure init_fun_interp_PRZ
end interface fun_interp_PRZ

!> interface for the mesh generation method
interface generate_st_mesh
  module procedure generate_random_mesh_1d
  module procedure generate_equidistant_mesh_1d
end interface generate_st_mesh

!> interface for the mesh generation method
interface generate_mesh_1d
  module procedure generate_random_mesh_1d
  module procedure generate_equidistant_mesh_1d
end interface generate_mesh_1d
contains
!> Procedures --------------------------------------------
!> generate random mesh in the s,t,phi coordinates
!> inputs:
!>   n_trials: (integer) number of mesh elements
!>   interval: (real8)(2) minimum and maximum mesh nodes
!>   rngs:     (type_rng)(n_threads) random number generators
!> outputs:
!>   rngs:     (type_rng)(n_threads) random number generators
!>   mesh:     (real8)(n_trials) mesh within intervals
subroutine generate_random_mesh_1d(n_trials,rngs,interval,mesh)
  use mod_rng, only: type_rng
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  class(type_rng),dimension(:),allocatable,intent(inout) :: rngs
  !> inputs:
  integer,intent(in)                                     :: n_trials
  real*8,dimension(2),intent(in)                         :: interval 
  !> outputs:
  real*8,dimension(n_trials),intent(out)                 :: mesh
  !> variables
  integer :: ii,n_trials_thread,thread_id,n_max_threads
  n_max_threads = 1;
  !$ n_max_threads = omp_get_max_threads()
  n_trials_thread = n_trials/n_max_threads
  !> generate mesh
  !$omp parallel do default(private) firstprivate(n_max_threads,&
  !$omp n_trials_thread,interval) shared(rngs,mesh)
  do ii=1,n_max_threads
    thread_id = 1;
    !$thread_id = omp_get_thread_num()
    call rngs(thread_id)%next(mesh((ii-1)*n_trials_thread+1:ii*n_trials_thread))
  enddo
  !$omp end parallel do
  if(n_trials.lt.n_max_threads*n_trials_thread) &
  call rngs(1)%next(mesh(n_trials_thread*n_max_threads+1:n_trials))
  mesh = interval(1)+(interval(2)-interval(1))*mesh
end subroutine generate_random_mesh_1d

!> generate equidistant mesh within interval
!> inputs:
!>   n_nodes:  (integer) number of mesh nodes
!>   interval: (real8)(2) first and last mesh nodes
!> outputs:
!>   mesh:  (real8)(n_nodes) equidistant mesh within interval
subroutine generate_equidistant_mesh_1d(n_nodes,interval,mesh)
  use constants, only: TWOPI
  implicit none
  !> inputs:
  integer,intent(in)             :: n_nodes
  real*8,dimension(2),intent(in) :: interval
  !> outputs:
  real*8,dimension(n_nodes),intent(out) :: mesh
  !> variables:
  integer :: ii,jj,kk
  real*8 :: d_node
  !> initialisation
  d_node   = (interval(2)-interval(1))/real(n_nodes-1,kind=8);
  !> generate equidistant points
  !$omp parallel do simd default(private) firstprivate(n_nodes,d_node,interval) &
  !$omp shared(mesh)
  do ii=1,n_nodes
    mesh(ii) = interval(1)+d_node*real(ii-1,kind=8)
  enddo
  !$omp end parallel do simd
end subroutine generate_equidistant_mesh_1d

!> find local and global minimum and maximum values for a jorek axisymmetric field
 !> inputs:
!>   field_id:            (integer) index of the jorek field
!>   n_mesh:              (integer) number of mesh element
!>   mesh:                (real8)(2,n_mesh) mesh element list (s,t)
!>   f_interp_PRZ_axisym: (fun_interp_PRZ_axisym) interpolation function for root finding
!> outputs:
!>   f_interp_PRZ_axisym: (fun_interp_PRZ_axisym) interpolation function for root finding
!>   min_list:            (real8)(n_elements) minimum of each element
!>   max_list:            (real8)(n_elements) maximum of each element
!>   minmax_global:       (real8)(2) global minimum and maximum of the field
!subroutine 

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
subroutine field_minmax(field_id,n_s,n_t,n_phi,mesh_s,mesh_t,&
mesh_phi,f_interp_PRZ,min_list,max_list,minmax_global)
  use mod_rootfinding,    only: newtons_method
  use mod_math_operators, only: compute_eigenvalues
  implicit none
  !> inputs-outputs:
  type(fun_interp_PRZ),intent(inout) :: f_interp_PRZ
  !> inputs:
  integer,intent(in)                 :: field_id,n_s,n_t,n_phi
  real*8,dimension(n_s),intent(in)   :: mesh_s
  real*8,dimension(n_t),intent(in)   :: mesh_t
  real*8,dimension(n_phi),intent(in) :: mesh_phi
  !> outputs:
  real*8,dimension(2),intent(out)       :: minmax_global
  real*8,dimension(f_interp_PRZ%element_list%n_elements),intent(out) :: min_list
  real*8,dimension(f_interp_PRZ%element_list%n_elements),intent(out) :: max_list
  !> variables:
  integer                 :: ii,jj,kk,pp,ierr,maxit,n_elements
  real*8                  :: tol
  real*8,dimension(3)     :: x_extrema,values
  real*8,dimension(3,3)   :: Jac
  complex*16,dimension(3) :: eigv
  !> initialisation
  tol=5.d-16; maxit=10000; n_elements=f_interp_PRZ%element_list%n_elements
  minmax_global = (/1.d21,-1.d21/); min_list = 1.d21; max_list = -1.d21;
  !> find minimum and maximum
  !$omp parallel do default(private) firstprivate(n_elements,n_s,n_t,n_phi,&
  !$omp maxit,tol,field_id) shared(f_interp_PRZ,mesh_s,mesh_t,mesh_phi) & 
  !$omp reduction(min:min_list) reduction(max:max_list) collapse(4)
  do ii=1,n_elements
    do jj=1,n_phi
      do kk=1,n_t
        do pp=1,n_phi
          !> find extrema
          call newtons_method(f_interp_PRZ,(/0.d0,0.d0,0.d0/),(/mesh_s(pp),&
          mesh_t(kk),mesh_phi(jj)/),x_extrema,2,(/ii,field_id/),ierr,tol,maxit)
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
    enddo
  enddo
  !$omp end parallel do
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

!> initialise the fun_interp_PRZ_axisym class
!> inputs:
!>   node_list:    (type_node_list) jorek node list
!>   element_list: (type_element_list) jorek element list
!> outputs:
!>   this: (fun_interp_PRZ_axisym) initialised fun_interp_PRZ_axisym class
function init_fun_interp_PRZ_axisym(node_list,element_list) result(this)
  use data_structure, only: type_node_list
  use data_structure, only: type_element_list
  implicit none
  !> inputs:
  type(type_node_list),intent(in)    :: node_list
  type(type_element_list),intent(in) :: element_list
  !> output:
  type(fun_interp_PRZ_axisym) :: this
  !> initialise class
  this%node_list = node_list; this%element_list = element_list;
end function init_fun_interp_PRZ_axisym

!> function used for computing the derivatives and the jacobian
!> of the derivatives of a jorek axisymmetric field (phi=0.d0)
!> inputs:
!>   this:         (fun_interp_PRZ_axisym) fun_interp_PRZ_axisym class
!>   x:            (real8)(2) (s,t) coordinates
!>   n_int_coords: (integer) must be 2
!>   int_coords:   (integer)(n_int_coords) integer coordinates:
!>                 1: element number
!>                 2: field number
!> outputs:
!>   this: (fun_interp_PRZ) fun_interp_PRZ class
!>   f:    (real8)(2) jorek field derivatives
!>   J:    (real8)(2,2) jorek field hessian
pure subroutine f_df_interp_PRZ_axisym(this,f,J,x,n_int_coords,int_coords)
  use mod_interp, only: interp_PRZ
  implicit none
  !> inputs-outputs:
  class(fun_interp_PRZ_axisym),intent(inout) :: this
  !> inputs:
  integer,intent(in)                         :: n_int_coords
  integer,dimension(n_int_coords),intent(in) :: int_coords
  real*8,dimension(2),intent(in)             :: x
  !> outputs:
  real*8,dimension(2),intent(out)   :: f
  real*8,dimension(2,2),intent(out) :: J
  !> variables:
  real*8 :: R,R_s,R_t,R_st,R_ss,R_tt
  real*8 :: Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
  real*8,dimension(1) :: P,P_s,P_t,P_phi
  real*8,dimension(1) :: P_st,P_ss,P_tt,P_sphi,P_tphi,P_phiphi

  !> interpolate the JOREK fields
  call interp_PRZ(this%node_list,this%element_list,int_coords(1),&
  (/int_coords(2)/),1,x(1),x(2),0.d0,P,P_s,P_t,P_phi,P_st,P_ss,P_tt,&
  P_sphi,P_tphi,P_phiphi,R,R_s,R_t,R_st,R_ss,R_tt,&
  Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
  !> fill up the values and the jacobian with the first
  !> and second order derivatives and store the value
  this%val = P(1); f = (/P_s,P_t/)
  J(:,1) = (/P_ss,P_st/)
  J(:,2) = (/P_st,P_tt/)
end subroutine f_df_interp_PRZ_axisym

!>--------------------------------------------------------
end module mod_fields_minmax
