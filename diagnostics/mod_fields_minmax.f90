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
subroutine field_minmax_monte_carlo(node_list,element_list,field_id,n_trials,&
phi_int,rng_type,rngs,minmax_list,minmax_global,my_id,n_tasks,ifail)
  use mod_settings,   only: n_vertex_max
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
  integer,intent(in)                                     :: field_id,n_trials
  integer,intent(in)                                     :: my_id,n_tasks
  real*8,dimension(2),intent(in)                         :: phi_int
  !> outputs:
  real*8,dimension(2),intent(out)                         :: minmax_global
  real*8,dimension(2,element_list%n_elements),intent(out) :: minmax_list
  !> variables:
  integer :: ii,jj,thread_id,n_vertex,n_trials_per_task
  real*8 :: R,Z
  real*8,dimension(1) :: val
  real*8,dimension(3) :: stphi_coords
  real*8,dimension(element_list%n_elements) :: min_list,max_list
  real*8,dimension(:,:),allocatable :: min_list_all,max_list_all
  !> initialisations
  if(allocated(rngs)) deallocate(rngs)
  call setup_shared_rngs(size(stphi_coords),rng_type,rngs)
  n_vertex = n_vertex_max; n_trials_per_task = n_trials/n_tasks
  if(my_id.eq.0) n_trials_per_task = n_trials - n_trials_per_task*(n_tasks-1)
  min_list = 1.d21; max_list = -1.d21;
  !$omp parallel default(private) firstprivate(n_trials_per_task,n_vertex,&
  !$omp phi_int,field_id) shared(node_list,element_list,rngs) &
  !$omp reduction(min:min_list) reduction(max:max_list)
  !> initialise the min_list and max_list using the element nodes
  !> apply the brute force method
  thread_id=1;
  !$ thread_id = omp_get_thread_num()+1
  !$omp do collapse(2)
  do ii=1,element_list%n_elements
    do jj=1,n_trials_per_task
      call rngs(thread_id)%next(stphi_coords)
      stphi_coords(1:2) = max(min(-1.d-1 + 1.2d0*stphi_coords(1:2),1.d0),0.d0)
      stphi_coords(3) = phi_int(1)+(phi_int(2)-phi_int(1))*stphi_coords(3)
      call interp_PRZ(node_list,element_list,ii,(/field_id/),1,&
      stphi_coords(1),stphi_coords(2),stphi_coords(3),val,R,Z)
      min_list(ii) = min(min_list(ii),val(1)); max_list(ii) = max(max_list(ii),val(1));     
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  !> compute the global minimum and maximum
  !> recover data from all mpi task and find the final minimum and maximum
  if(my_id.eq.0) then
    allocate(min_list_all(element_list%n_elements,n_tasks)); 
    allocate(max_list_all(element_list%n_elements,n_tasks))
    call MPI_Gather(min_list,element_list%n_elements,MPI_DOUBLE,min_list_all,&
    element_list%n_elements,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
    call MPI_Gather(max_list,element_list%n_elements,MPI_DOUBLE,max_list_all,&
    element_list%n_elements,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  else
    call MPI_Gather(min_list,element_list%n_elements,MPI_DOUBLE,min_list_all,&
    0,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
    call MPI_Gather(max_list,element_list%n_elements,MPI_DOUBLE,max_list_all,&
    0,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  endif
  if(my_id.eq.0) then
    !$omp parallel do default(private) shared(element_list,min_list_all,&
    !$omp max_list_all,minmax_list)
    do ii=1,element_list%n_elements
      minmax_list(:,ii) = (/minval(min_list_all(ii,:)),maxval(max_list_all(ii,:))/)
    enddo
    !$omp end parallel do
    deallocate(max_list_all); deallocate(min_list_all);
  endif
  call MPI_BCast(minmax_list,2*element_list%n_elements,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  minmax_global(1)=minval(minmax_list(1,:)); minmax_global(2)=maxval(minmax_list(2,:))  
end subroutine field_minmax_monte_carlo

!>--------------------------------------------------------
end module mod_fields_minmax
