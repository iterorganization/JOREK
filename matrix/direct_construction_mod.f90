module direct_construction_mod

implicit none
  logical, private                   :: matrix_structure_initialized = .false.  

contains  

  !> Constructing harmonic matrix directly from the elementary matrix 
  subroutine direct_construction_harmonic(my_id, my_id_n, m_cpu, n_cpu, MPI_COMM_N,  MPI_COMM_MASTER, my_id_master, & 
    node_list, element_list, bnd_elm_list, bnd_node_list, xpoint2, xcase2, freeboundary, direct_construction)

  use data_structure 
  use global_distributed_matrix
  use mod_global_matrix_structure
  use equil_info
  use construct_matrix_mod, only : construct_matrix 
  use mpi_mod
  use vacuum
  use mod_integer_types
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in)    :: node_list
  type(type_bnd_node_list),     intent(inout) :: bnd_node_list
  type (type_element_list),     intent(in)    :: element_list
  type (type_bnd_element_list), intent(in)    :: bnd_elm_list
  integer, intent(in) :: my_id, my_id_n, n_cpu, m_cpu, MPI_COMM_N, MPI_COMM_MASTER, my_id_master, xcase2
  logical, intent(in) :: direct_construction, xpoint2, freeboundary
  
  ! --- Local variables
  integer(kind=int_all)              :: Int1=1
  integer,               allocatable :: index_min_harm(:), index_max_harm(:)
  integer,               allocatable :: local_elms_harm(:)
  integer                            :: n_local_elms_harm
  integer(kind=int_all)              :: ndof 
  integer                            :: i_tor_min, i_tor_max 
  integer                            :: i, ierr
    
  ! --- Memory allocation 
  if (allocated(local_elms_harm)) call tr_deallocate(local_elms_harm,"local_elms_harm",CAT_DMATRIX) 
  if (allocated(index_min_harm))  call tr_deallocate(index_min_harm,"index_min_harm",CAT_DMATRIX) 
  if (allocated(index_max_harm))  call tr_deallocate(index_max_harm,"index_max_harm",CAT_DMATRIX) 
 
  call tr_allocate(local_elms_harm,1,element_list%n_elements,"local_elms_harm",CAT_FEM)
  call tr_allocate(index_min_harm,1,n_cpu,"index_min_harm",CAT_FEM)
  call tr_allocate(index_max_harm,1,n_cpu,"index_max_harm",CAT_FEM)

  if(my_id .lt. m_cpu)  then
    i_tor_min = 1
    i_tor_max = 1
  else
    i_tor_min = 2*(my_id - MOD(my_id, m_cpu))/m_cpu
    i_tor_max = i_tor_min + 1
  endif
 
  call distribute_nodes_elements(my_id,m_cpu,n_cpu,node_list,element_list, direct_construction, & 
    local_elms_harm, n_local_elms_harm, ndof, index_min_harm,index_max_harm)

  if ( .not. matrix_structure_initialized ) then
    matrix_structure_initialized = .true.

    call global_matrix_structure(my_id,my_id_n,node_List,element_list,bnd_elm_list, freeboundary, &
      local_elms_harm,n_local_elms_harm,index_min_harm(my_id+1),                                  & 
      index_max_harm(my_id+1), ijA_index_harm, ijA_size_harm,                                     &
      irn_jcn_harm, irn_harm, jcn_harm, i_tor_min, i_tor_max,                           &                         
      n_harm, nz_harm, ndof_harm, n_matrix_block_size_harm)

    if ( freeboundary .and. ( sr%n_tor /= 0 ) ) then 
      call global_matrix_structure_vacuum(node_list, bnd_node_list, index_min_harm(my_id+1), index_max_harm(my_id+1), & 
        i_tor_min, i_tor_max, irn_harm, jcn_harm, n_matrix_block_size_harm, ijA_index_harm, ijA_size_harm, irn_jcn_harm) 
    endif
  endif
 

  call construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master,                &
    local_elms_harm, n_local_elms_harm, index_min_harm(my_id+1), index_max_harm(my_id+1), xpoint2,&
    xcase2, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%psi_bnd, ES%R_xpoint, ES%Z_xpoint,              &
    ES%psi_xpoint, i_tor_min, i_tor_max, n_harm, nz_harm, ndof_harm, n_matrix_block_size_harm, A_harm,  &
    rhs_harm, irn_harm, jcn_harm, ijA_index_harm, ijA_size_harm, irn_jcn_harm,     &
    direct_construction)

  end subroutine direct_construction_harmonic

end module direct_construction_mod
