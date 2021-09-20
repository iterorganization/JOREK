module preconditioner_module
  use phys_module, only: modes_per_family, mode_families_modes, autodistribute_modes, n_mode_families, weights_per_family, &
                         autodistribute_ranks, ranks_per_family, centralize_harm_mat, use_strumpack
  use mod_integer_types

  implicit none

  integer :: my_family_id
  integer, dimension(:), allocatable :: my_mode_set !< Mode number in local mode family used for preconditioner
  integer(kind=int_all), dimension(:), allocatable :: my_row_index !< Row indices of local mode family in global RHS  - can be replaced by logical
!  real, dimension(:), allocatable :: my_row_factor !< Multiplying factor of local mode family in global RHS  - can be replaced by logical
  real :: my_row_factor !< Multiplying factor of local mode family in global RHS  - can be replaced by logical
  integer :: my_mode_set_n !< number of modes in local mode family
  integer, allocatable :: mode_families_ranks(:,:), rank_range(:)

  private
  public create_communicators, distribute_modes, check_preconditioner_consistency, &
         my_family_id, my_row_index, my_mode_set, my_mode_set_n, &
         map_row_index, my_row_factor, rank_range, mode_families_ranks

  contains

  !> Distribute MPI ranks among mode families
  !! i_tor(n_cpu) provides family ID for each rank within MPI_COMM_WORLD
  subroutine distribute_ranks(n_cpu,i_tor)

    implicit none

    integer, intent(in) :: n_cpu
    integer, intent(out) :: i_tor(n_cpu)
    integer :: mcpu, r, i, j

    mcpu = n_cpu/n_mode_families
    r = mod(n_cpu,n_mode_families)

    allocate(rank_range(n_mode_families+1))
    if (allocated(mode_families_ranks)) deallocate(mode_families_ranks)
    allocate(mode_families_ranks(n_mode_families,n_cpu))
    mode_families_ranks = -1

    i_tor = 0

    if (autodistribute_ranks) then
      do i=1, n_mode_families
        ranks_per_family(i) = mcpu
        if ((r.gt.0).and.(i.le.r))  ranks_per_family(i) = ranks_per_family(i) + 1 ! add extra rank if avaiable
      enddo
    endif

    rank_range(1) = 1
    do i = 2, n_mode_families+1
      rank_range(i) = rank_range(i-1) + ranks_per_family(i-1)
    enddo

    ! check for consistency
    r = 0
    do i=2,n_mode_families+1
      r = r + rank_range(i) - rank_range(i-1)
    enddo
    if (r.ne.n_cpu) then
      write(*,*) "Error in distribution of ranks"
      call exit(0)
    endif

    do i=1,n_cpu
      do j=2,n_mode_families+1
        if ((i.ge.rank_range(j-1)).and.(i.lt.rank_range(j))) then
          i_tor(i) = j-1
          exit
        endif
      enddo
    enddo

    do j=1,n_mode_families
      r = 0
      do i = 1,n_cpu
        if (i_tor(i).eq.j) then
          r = r + 1
          mode_families_ranks(j,r) = i - 1
        endif
      enddo
    enddo

  end subroutine distribute_ranks

  !> Distribute toroidal modes among mode families
  subroutine distribute_modes

    use mod_parameters, only : n_tor
    implicit none

    integer :: i

    if (autodistribute_modes) then
      modes_per_family(1) = 1
      mode_families_modes(1,1) = 1
      if (n_mode_families.gt.1) then
        do i = 2, n_mode_families
          modes_per_family(i) = 2
          mode_families_modes(i,1) =  (i - 1)*2
          mode_families_modes(i,2) =  (i - 1)*2 + 1
        enddo
      endif
    endif

    if (allocated(my_mode_set)) deallocate(my_mode_set)
    my_mode_set_n = modes_per_family(my_family_id)
    allocate(my_mode_set(my_mode_set_n))
    my_mode_set(1:my_mode_set_n) = mode_families_modes(my_family_id,1:my_mode_set_n)

    write(*,*) "my_family_id:", my_family_id, "my_mode_set:", my_mode_set
 
  end subroutine distribute_modes

  !> Set up MPI communicators for mode families and corresponding masters
  subroutine create_communicators(my_id_n, n_cpu_n, MPI_COMM_N, my_id_master, n_masters, MPI_COMM_MASTER, MPI_COMM_TRANS)

    use mpi
    implicit none

    integer, intent(out) :: my_id_n, my_id_master, n_cpu_n, n_masters, MPI_COMM_N, MPI_COMM_MASTER, MPI_COMM_TRANS

    integer, allocatable :: i_tor(:)
    integer :: i, my_id, n_cpu, ierr, MPI_GROUP_WORLD, MPI_GROUP_MASTER
    integer, allocatable :: ranks_tmp(:)

    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)

    if (n_cpu.lt.n_mode_families) then
      write(*,*) "Error: number of ranks must be >= n_mode_families"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif

    allocate(i_tor(n_cpu))
    call distribute_ranks(n_cpu,i_tor)
    if (my_id.eq.0) write(*,*) "ranks_per_family", ranks_per_family(1:n_mode_families)

    my_family_id = i_tor(my_id+1)

    call MPI_COMM_SPLIT(MPI_COMM_WORLD,i_tor(my_id+1),my_id,MPI_COMM_N,ierr)
    if (ierr.ne.0) then
      write(*,*) "Error in creating MPI_COMM_N"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif
    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)
    call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)

    i_tor(my_id+1) = my_id_n
    call MPI_Allreduce(MPI_IN_PLACE,i_tor,n_cpu,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,i_tor(my_id+1),my_id,MPI_COMM_TRANS,ierr)

    n_masters = n_mode_families
    allocate(ranks_tmp(n_masters)); ranks_tmp=0;
    if (my_id_n.eq.0) ranks_tmp(my_family_id) = my_id
    call MPI_AllReduce(MPI_IN_PLACE,ranks_tmp,n_masters,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
    call MPI_GROUP_INCL(MPI_GROUP_WORLD,n_masters,ranks_tmp,MPI_GROUP_MASTER,ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP_MASTER,MPI_COMM_MASTER,ierr)
    deallocate(ranks_tmp)

    if (my_id_n .eq. 0) then
     call MPI_COMM_RANK(MPI_COMM_MASTER, my_id_master, ierr)
    endif

    if ((my_id.eq.0).and.(my_id_n.ne.0)) then
      write(*,*) "Error in creating communicators: my_id==0 must have my_id_n==0"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif

    deallocate(i_tor)

    if (my_id.eq.0) then
      do i=1, n_mode_families
        write(*,*) "mode_families_ranks", i, mode_families_ranks(i,1:ranks_per_family(i))
      enddo
    endif

    return

  end subroutine create_communicators

  !> Determine mapping from local to globar row index for the RHS
  subroutine map_row_index(ndof)

    use mod_parameters, only: n_tor
    use tr_module
    use mod_integer_types

    implicit none

    integer(kind=int_all), intent(in) :: ndof
    integer(kind=int_all) :: i, ndof_family
    integer(kind=int_all), parameter   :: Int1=1
    integer :: im

    ndof_family = my_mode_set_n*ndof/n_tor

    if (allocated(my_row_index)) call tr_deallocate(my_row_index,"my_row_index",CAT_DMATRIX)
    call tr_allocate(my_row_index,Int1,ndof_family,"my_row_index",CAT_DMATRIX)
    !if (allocated(my_row_factor)) call tr_deallocate(my_row_factor,"my_row_index",CAT_DMATRIX)
    !call tr_allocate(my_row_factor,Int1,ndof_family,"my_row_factor",CAT_DMATRIX)
    my_row_factor = weights_per_family(my_family_id)

    do i = 0, ndof/n_tor - 1
      do im=1,my_mode_set_n
        my_row_index(im+i*my_mode_set_n) =  my_mode_set(im) + i*n_tor
      enddo
    enddo

    if (autodistribute_modes) my_row_factor = 1.0

    write(*,*) "my_family_id my_row_factor", my_family_id, weights_per_family(my_family_id)

    return
  end subroutine map_row_index

  !> Check mode families consistency
  subroutine check_preconditioner_consistency

    use mpi
    use mod_parameters, only: n_tor

    implicit none

    integer :: i, j, my_id, n_cpu, ierr
    logical :: cck

    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)

    if (n_mode_families.gt.n_cpu) then
      if (my_id.eq.0) write(*,*) "Error: number of cpu must be >= number of mode families", n_mode_families
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif

    if (.not.autodistribute_modes) then
      ! check if mode families are specified
      if (n_mode_families.le.0) then
        if (my_id.eq.0) write(*,*) "Error: number of mode families is not specified or <=0"
        call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
      endif

      ! check if mode families consist of correct modes
      do i = 1, n_mode_families
        do j = 1, modes_per_family(i)
          if ((mode_families_modes(i,j).lt.1).or.(mode_families_modes(i,j).gt.n_tor)) then
            if (my_id.eq.0) write(*,*) "Error: incorrect specification of mode_families_modes", i
            call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
          endif
        enddo
      enddo

    endif

    if (.not.autodistribute_ranks) then
      ! check if numbers of ranks are provided
      do i = 1, n_mode_families
        if (ranks_per_family(i).le.0) then
          if (my_id.eq.0) write(*,*) "Error: ranks_per_family is not correctly specified", i
          call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
        endif
      enddo
      ! check if sum of all ranks per family is equal to the total number of ranks
      if (sum(ranks_per_family(1:n_mode_families)).ne.n_cpu) then
        if (my_id.eq.0) write(*,*) "Error: sum of all ranks_per_family must be equal to total number of MPI ranks"
        call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
      endif
    endif

    if ((.not.use_strumpack).and.(.not.centralize_harm_mat)) then
      if (my_id.eq.0) write(*,*) "Warning: PC matrix distribution is only supported by STRUMPACK"
      centralize_harm_mat = .true.
    elseif ((use_strumpack).and.(centralize_harm_mat)) then
      if (my_id.eq.0) write(*,*) "Warning: centralization of PC matrix is not advized when using STRUMPACK"
    endif


  end subroutine check_preconditioner_consistency

end module preconditioner_module
