!> Module allowing to exchange indices (node%index) of the finite element grid in order to optimize
!! the parallelization of the free boundary part of the code. The exchanging of the indices makes
!! sure that all MPI tasks contribute to the calculation of the boundary integral, i.e., feel
!! responsible for some of the boundary elements.
module mod_exchange_indices

implicit none

logical, parameter         :: DEBUG_OUTPUT      = .false. !< Hard-coded parameter for debug output
logical, save              :: initialized       = .false. !< Has the module been initialized?
logical, save              :: indices_exchanged = .false. !< Have the indices been exchanged w.r.t.
                                                          !! their normal order?
integer, save              :: n_nodes = -99               !< Stored value of node_list%n_nodes
integer, save              :: len_exchange                !< How many entries in the exchange table?
integer, save, allocatable :: exchange_table(:,:)         !< Table with the indices that should be
                                                          !! exchanged
integer, allocatable :: index_min(:), index_max(:)        !< Min/max indices for each MPI rank



private
public exchange_indices



contains



!> Is a specific "my_id" responsible for a given node index?
logical function is_responsible(index, irank)
  
  integer, intent(in) :: index, irank
  
  is_responsible = ( index >= index_min(irank) ) .and. ( index <= index_max(irank) )
  
end function is_responsible



!> Add another pair of exchange rules to the table
subroutine add_exchange_rules(ind1, ind2, num)
  
  integer, intent(in)    :: ind1, ind2
  integer, intent(inout) :: num
  
  integer :: i
  
  do i = 1, num-1
    if ( exchange_table(i,1) == ind2 ) return
  end do
  
  exchange_table(num,:) = (/ind1, ind2/)
  num = num + 1
  exchange_table(num,:) = (/ind2, ind1/)
  num = num + 1
  
end subroutine add_exchange_rules



!> Exchange some indices of grid nodes in order to parallelize the vacuum boundary integral.
subroutine initialize(node_list, my_id, n_cpu)
  
  use data_structure
  
  ! --- Routine parameters
  type(type_node_list), intent(inout) :: node_list
  integer,              intent(in)    :: my_id
  integer,              intent(in)    :: n_cpu

  ! --- Local variables
  integer, allocatable :: first_index_usable(:), mm(:)
  integer :: i, j, k, l, ind_max, n_bnd, ind_bnd, ind1, ind2
  
  if ( DEBUG_OUTPUT ) write(*,*) 'Initializing module mod_exchange_indices'
  
  ! --- Determine maximum index in the grid and number of boundary nodes
  ind_max = -1
  n_bnd   = 0
  do i = 1, node_list%n_nodes
    ind_max = max(ind_max, maxval(node_list%node(i)%index(:)))
    if ( node_list%node(i)%boundary > 0 ) n_bnd = n_bnd + 1
  end do
  if ( DEBUG_OUTPUT .and. (my_id == 0) ) then
    write(*,*) 'n_cpu   =', n_cpu
    write(*,*) 'ind_max =', ind_max
    write(*,*) 'n_bnd   =', n_bnd
  end if
  
  ! --- Determine the index_min and index_max locally in the module
  !     (as distribute_nodes_elements will be called only later on)
  if ( allocated(index_min) ) deallocate( index_min )
  if ( allocated(index_max) ) deallocate( index_max )
  allocate( index_min(n_cpu), index_max(n_cpu) )
  index_min(1) = 1
  do i = 1, n_cpu
    index_max(i) = (i * ind_max) / n_cpu
  enddo
  do i=2,n_cpu
    index_min(i) = index_max(i-1) + 1
  enddo
  if ( DEBUG_OUTPUT ) then
    write(*,'(a,99i7)') ' index_min =', index_min(:)
    write(*,'(a,99i7)') ' index_max =', index_max(:)
  end if
  
  ! --- Find out which grid nodes are usable for "exchanging indices"; store the number of the
  !     first index for each MPI rank that can be used; skip the grid center nodes as these might
  !     not have four independent grid indices
  if ( allocated(exchange_table) )     deallocate( exchange_table )
  if ( allocated(first_index_usable) ) deallocate( first_index_usable )
  if ( allocated(mm) )                 deallocate( mm )
  allocate(exchange_table(n_bnd*8,2))
  allocate(first_index_usable(n_cpu))
  allocate(mm(n_cpu))
  first_index_usable(:) = -1
  ii: do i = 1, n_cpu
    do j = 1, node_list%n_nodes
      if ( is_responsible(node_list%node(j)%index(1), i) .and. (node_list%node(j)%index(1)>1) ) then
        first_index_usable(i) = node_list%node(j)%index(1)
        cycle ii
      end if
    end do
  end do ii
  
  if ( minval(first_index_usable) < 1 ) then
    write(*,*) my_id, 'ERROR: first_index_usable < 1'
    write(*,*) my_id, first_index_usable(:)
    stop
  end if
  
  ! --- Prepare a "table" of indices to be exchanged
  ind_bnd = 0
  j       = 1
  mm(:)   = first_index_usable(:)
  if ( DEBUG_OUTPUT ) write(*,*) 'mm before:', mm(:)
  do i = 1, node_list%n_nodes
    if ( node_list%node(i)%boundary > 0 ) then
      k = (real(j)/real(8*n_bnd))*n_cpu + 1 ! with which MPI rank to we want to exchange this index?
      if ( k == n_cpu ) cycle ! the last MPI rank doesn't need to exchange with itself
      do l = 1, 4 ! the four dofs of one node
        ind1    = node_list%node(i)%index(l) ! exchange this index
        ind2    = mm(k) + l - 1              ! with this one for which MPI rank k is responsible
        if ( is_responsible(ind1,k) ) cycle ! need not exchange as MPI rank is already responsible
        call add_exchange_rules(ind1, ind2, j)
        if ( DEBUG_OUTPUT ) write(*,*) 'LIST: ', ind1, '<->', ind2
      end do
      mm(k) = mm(k) + 4
    end if
  end do
  len_exchange = j - 1
  if ( DEBUG_OUTPUT .and. (my_id == 0) ) write (*,*) 'len_exchange ', len_exchange
  
  deallocate(first_index_usable, mm, index_min, index_max)
  
  initialized = .true.
  
end subroutine initialize


  
!> Exchange some indices of grid nodes in order to parallelize the vacuum boundary integral.
subroutine exchange_indices(node_list, my_id, n_cpu, back)
  
  use data_structure
    
  ! --- Routine parameters
  type(type_node_list), intent(inout) :: node_list
  integer,              intent(in)    :: my_id
  integer,              intent(in)    :: n_cpu
  logical,              intent(in)    :: back   !< Change indices back (used to check alternating
                                                !! order of forth and back exchanges only)
  
  ! --- Local variables
  integer :: i, j, k, l
  
  ! --- Needs re-initialization due to grid change?
  if ( initialized .and. (n_nodes /= node_list%n_nodes) ) then
    write(*,*) 'The grid seems to have changed. Will re-initialize mod_exchange_indices.'
    initialized = .false.
  end if
  n_nodes = node_list%n_nodes
  
  ! --- A few checks and warnings
  if ( my_id == 0 ) write(*,*) 'EXCHANGE_INDICES to enhance free boundary simulation performance'
  if ( DEBUG_OUTPUT ) write(*,*) my_id, n_cpu, back
  
  if ( n_cpu == 1 ) then
    write(*,*) my_id, 'Remark: Exchange_indices is skipped as you are running with a single MPI task.'
    return
  else if ( indices_exchanged .neqv. back ) then
    write(*,*) my_id, 'ERROR: Somewhere in the code you call exchange_indices too often or not often enough.'
    stop
  else if ( indices_exchanged .and. ( .not. initialized ) ) then
    write(*,*) my_id, 'ERROR: Internal bug! indices_exchanged=.t. and initialized=.f. should not happen'
    stop
  end if
  
  ! --- Initialize when called the first time in a run (or after a grid change)
  if ( .not. initialized ) call initialize(node_list, my_id, n_cpu)
  
  ! --- Exchange the indices
  l = 0
  do i = 1, node_list%n_nodes
    do j = 1, 4
      do k = 1, len_exchange
        if ( node_list%node(i)%index(j) == exchange_table(k,1) ) then
          if ( DEBUG_OUTPUT ) write(*,'(a,i7,a,i7,a,i7,a,i3,a)') 'ex:', node_list%node(i)%index(j), '->', exchange_table(k,2), ' (node', i, ', dof', j, ')'
          node_list%node(i)%index(j) = exchange_table(k,2)
          l = l + 1
          exit
        end if
      end do
    end do
  end do
  if ( len_exchange /= l ) write(*,*) my_id, 'WARNING: The actual number of exchanged indices does not match the expected value.', l, len_exchange
  if ( DEBUG_OUTPUT .and. (my_id == 0) ) write (*,*) 'num exchanged ', l
  
  indices_exchanged = .not. indices_exchanged ! switch the state
  
end subroutine exchange_indices



end module mod_exchange_indices
