!> mod_field_minmax contains variables and procedures for 
!> for finding the local and global maximum and minimum
!> values of node list fields
module mod_fields_minmax
implicit none
private
public :: field_minmax_monte_carlo

!> Variables and datatypes -------------------------------
!> Interfaces---------------------------------------------
contains
!> try to find  local minimum and maximum via brute force method
!> inputs:
!>   node_list:    (type_node_list) jorek mesh node list
!>   element_list: (type_element_list) jorek mesh element list
!>   n_fields:     (integer) number of fluid fierlds to treat
!>   fields_ids:   (integer)(n_fields) node indices of each field
!>   n_trials:     (integer) number of trials of the testing method
!>   phi_int:      (real8)(2) toroidal angle interval
!>   rng_type:     (type_rng) type of the RNG to use
!>   rngs:         (type_rng)(:)(allocatable) array of the RNGs to use
!>   my_id:        (integer) id of the MPI task
!>   n_tasks:      (integer) total number of MPI tasks
!>   ifail:        (integer) MPI error
!> outputs:
!>   minmax_list:   (real8)(n_fields,2,n_elements) estimated minimum
!>                  and maximum for all fields for each element
!>   minmax_global: (real8)(n_fields,2) estimated global minimum
!>                  and maximum for all fields
!>   ifail:        (integer) MPI error
subroutine field_minmax_monte_carlo(node_list,element_list,n_fields,field_ids,&
n_trials,phi_int,rng_type,rngs,minmax_list,minmax_global,my_id,n_tasks,ifail)
  use data_structure, only: type_node_list,type_element_list
  use mod_interp,     only: interp_PRZ
  use mod_rng,        only: setup_shared_rngs
  use mod_rng,        only: type_rng
  use mpi
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  type(type_node_list),intent(in)                        :: node_list
  type(type_element_list),intent(in)                     :: element_list
  class(type_rng),dimension(:),allocatable,intent(inout) :: rngs
  class(type_rng),intent(in) :: rng_type
  integer,intent(in)                                     :: n_fields,n_trials
  integer,intent(in)                                     :: my_id,n_tasks
  integer,dimension(n_fields),intent(in)                 :: field_ids
  real*8,dimension(2),intent(in)                         :: phi_int
  !> outputs:
  real*8,dimension(n_fields,2),intent(out)               :: minmax_global
  real*8,dimension(n_fields,2,element_list%n_elements),intent(out) :: minmax_list
  !> variables:
  integer :: ii,jj,thread_id,n_trials_per_task
  real*8 :: R,Z
  real*8,dimension(3) :: stphi_coords
  real*8,dimension(n_fields) :: vals
  real*8,dimension(n_fields,element_list%n_elements) :: min_list,max_list
  !> initialisations
  if(allocated(rngs)) deallocate(rngs); thread_id=1;
  call setup_shared_rngs(size(stphi_coords),rng_type,rngs)
  n_trials_per_task = n_trials/n_tasks
  if(my_id.eq.0) n_trials_per_task = n_trials - n_trials_per_task*(n_tasks-1)
  !$omp parallel default(private) firstprivate(n_trials_per_task,phi_int,&
  !$omp n_fields,thread_id) shared(node_list,element_list,field_ids,rngs,&
  !$omp min_list,max_list)
  !$ thread_id = omp_get_thread_num()+1
  !$omp do collapse(2) reduction(min:min_list) reduction(max:max_list)
  do jj=1,n_trials_per_task
    do ii=1,element_list%n_elements
      call rngs(thread_id)%next(stphi_coords)
      stphi_coords(1:2) = max(min(-1.d-1 + 1.2d0*stphi_coords(1:2),1.d0),0.d0)
      stphi_coords(3) = phi_int(1)+(phi_int(2)-phi_int(1))*stphi_coords(3)
      call interp_PRZ(node_list,element_list,ii,field_ids,n_fields,&
      stphi_coords(1),stphi_coords(2),stphi_coords(3),vals,R,Z)
      min_list(:,ii) = min(min_list(:,ii),vals) 
      max_list(:,ii) = max(max_list(:,ii),vals)
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  !> compute the global minimum and maximum
  if(my_id.eq.0) then
    call MPI_Reduce(MPI_IN_PLACE,min_list,n_fields*2,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ifail)
    call MPI_Reduce(MPI_IN_PLACE,max_list,n_fields*2,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ifail)
    minmax_list(:,1,:) = min_list; minmax_list(:,2,:) = max_list;
  else
    call MPI_Reduce(min_list,min_list,n_fields*2,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ifail)
    call MPI_Reduce(max_list,max_list,n_fields*2,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ifail)
  endif 
  call MPI_BCast(minmax_list,n_fields*2*element_list%n_elements,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  !> estimate the global maximum and minimum
  minmax_global(:,1) = minval(minmax_list(:,1,:),dim=2)
  minmax_global(:,2) = maxval(minmax_list(:,2,:),dim=2)
end subroutine field_minmax_monte_carlo

!>--------------------------------------------------------
end module mod_fields_minmax
