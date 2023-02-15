module mod_preconditioner

  private
  public initialize_preconditioner, reset_preconditioner, update_pc_rhs, gather_solution

  contains

!> Initialize preconditioner structure (PC)
!! call subroutines for setting mode families, creating communicators, distributing tasks
  subroutine initialize_preconditioner(pc,comm_glob)
    use phys_module, only: autodistribute_modes, n_mode_families, autodistribute_ranks, centralize_harm_mat
    use mod_parameters, only: n_tor
    use data_structure, only: type_PRECOND
    use mpi_mod
    implicit none

    type(type_PRECOND) :: pc
    integer            :: comm_glob, my_id, n_cpu, ierr
    integer            :: i
    character(len=256) :: s

    pc%comm = comm_glob

    call MPI_COMM_RANK(pc%comm, my_id, ierr)
    call MPI_COMM_SIZE(pc%comm, n_cpu, ierr)

    pc%my_id = my_id
    pc%n_cpu = n_cpu
    if (pc%my_id.eq.0) write(*,*) "Initializing preconditioner"

    pc%autodistribute_ranks = autodistribute_ranks
    pc%autodistribute_modes = autodistribute_modes
    !pc%mat%row_distributed  = .not.centralize_harm_mat

    if (pc%autodistribute_modes) then
      pc%n_mode_families = (n_tor + 1)/2
    else
      pc%n_mode_families = n_mode_families
    endif

    call distribute_ranks(n_cpu, pc)

    call create_communicators(pc)
    pc%mat%comm = pc%MPI_COMM_N ! communicator for PC matrix distribution

    call distribute_modes(pc)

    pc%initialized = .true.

    if (my_id.eq.0) then
      do i=1, pc%n_mode_families
        write(s,'(A17,i4,A12)') " mode_family_id: ", i, " MPI ranks: "
        write(*,*) trim(s), pc%mode_families_ranks(i,1:pc%ranks_per_family(i))
      enddo
      do i=1, pc%n_mode_families
        write(s,'(A17,i4,A9,f6.2)') " mode_family_id: ", i, " weight: ", pc%row_factor
        write(*,*) trim(s), " modes:", pc%mode_families_modes(i,1:pc%modes_per_family(i))
      enddo
    endif

    return

  end subroutine initialize_preconditioner

!> Set up MPI communicators for mode families and corresponding masters
  subroutine create_communicators(pc)
    use data_structure, only: type_PRECOND
    use mpi_mod
    implicit none

    type(type_PRECOND) :: pc

    integer, allocatable :: i_tor(:), ranks_tmp(:)
    integer :: my_id, n_cpu, ierr

    my_id = pc%my_id
    n_cpu = pc%n_cpu

    call MPI_COMM_SPLIT(pc%comm, pc%family_id, my_id, pc%MPI_COMM_N, ierr)
    if (ierr.ne.0) then
      write(*,*) "Error in creating MPI_COMM_N"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif
    call MPI_COMM_RANK(pc%MPI_COMM_N, pc%my_id_n, ierr)
    call MPI_COMM_SIZE(pc%MPI_COMM_N, pc%n_cpu_n, ierr)

    allocate(i_tor(n_cpu)); i_tor = 0
    i_tor(my_id + 1) = pc%my_id_n
    call MPI_Allreduce(MPI_IN_PLACE,i_tor,n_cpu,MPI_INT,MPI_SUM,pc%comm,ierr)
    call MPI_COMM_SPLIT(pc%comm,i_tor(my_id+1),my_id,pc%MPI_COMM_TRANS,ierr)

    pc%n_masters = pc%n_mode_families
    allocate(ranks_tmp(pc%n_masters)); ranks_tmp=0;

    if (pc%my_id_n.eq.0) ranks_tmp(pc%family_id) = my_id
    call MPI_AllReduce(MPI_IN_PLACE,ranks_tmp,pc%n_masters,MPI_INT,MPI_SUM,pc%comm,ierr)
    call MPI_COMM_GROUP(pc%comm,pc%MPI_GROUP_WORLD,ierr)
    call MPI_GROUP_INCL(pc%MPI_GROUP_WORLD,pc%n_masters,ranks_tmp,pc%MPI_GROUP_MASTER,ierr)
    call MPI_COMM_CREATE(pc%comm,pc%MPI_GROUP_MASTER,pc%MPI_COMM_MASTER,ierr)

    if (pc%my_id_n .eq. 0) then
     call MPI_COMM_RANK(pc%MPI_COMM_MASTER, pc%my_id_master, ierr)
    endif

    if ((my_id.eq.0).and.(pc%my_id_n.ne.0)) then
      write(*,*) "Error in creating communicators: my_id==0 must have my_id_n==0"
      call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
    endif

    deallocate(i_tor, ranks_tmp)

    return

  end subroutine create_communicators

  !> Distribute toroidal modes among mode families
  subroutine distribute_modes(pc)
    use data_structure, only: type_PRECOND
    use phys_module, only: modes_per_family, mode_families_modes, weights_per_family
    implicit none

    type(type_PRECOND) :: pc
    integer            :: i, j, n_fam_max

    allocate(pc%modes_per_family(pc%n_mode_families))

    if (pc%autodistribute_modes) then
      pc%row_factor = 1.0
      pc%modes_per_family(1) = 1
      if (pc%n_mode_families>1) pc%modes_per_family(2:pc%n_mode_families) = 2
    else
      do i = 1, pc%n_mode_families
        pc%row_factor = weights_per_family(i)
        pc%modes_per_family(i) = modes_per_family(i)
      enddo
    endif

    n_fam_max = 1
    do i = 1, pc%n_mode_families
      n_fam_max = max(n_fam_max,pc%modes_per_family(i))
    enddo

    allocate(pc%mode_families_modes(pc%n_mode_families,n_fam_max))
    pc%mode_families_modes(:,:) = -1

    if (pc%autodistribute_modes) then
      pc%mode_families_modes(1,1) = 1
      if (pc%n_mode_families.gt.1) then
        do i = 2, pc%n_mode_families
          pc%mode_families_modes(i,1) =  (i - 1)*2
          pc%mode_families_modes(i,2) =  (i - 1)*2 + 1
        enddo
      endif
    else
      do i = 1, pc%n_mode_families
        do j = 1, modes_per_family(i)
          pc%mode_families_modes(i,j) = mode_families_modes(i,j)
        enddo
      enddo
    endif

    pc%mode_set_n = pc%modes_per_family(pc%family_id)
    allocate(pc%mode_set(pc%mode_set_n))
    pc%mode_set(1:pc%mode_set_n) = pc%mode_families_modes(pc%family_id,1:pc%mode_set_n)

  end subroutine distribute_modes

  !> Distribute MPI ranks among mode families
  subroutine distribute_ranks(n_cpu,pc)
    use data_structure, only: type_PRECOND
    use phys_module, only: ranks_per_family
    implicit none

    type(type_PRECOND) :: pc
    integer, intent(in)  :: n_cpu
    integer :: mcpu, r, i, j
    integer, dimension(:), pointer :: rank_id => Null()

    allocate(pc%rank_range(pc%n_mode_families + 1))
    allocate(pc%ranks_per_family(pc%n_mode_families))
    allocate(pc%mode_families_ranks(pc%n_mode_families,n_cpu))
    allocate(rank_id(n_cpu))

    do i = 1, pc%n_mode_families
      do j = 1, n_cpu
        pc%mode_families_ranks(i,j) = -1
      enddo
    enddo

    mcpu = n_cpu/pc%n_mode_families
    r = mod(n_cpu,pc%n_mode_families)

    if (pc%autodistribute_ranks) then
      do i=1, pc%n_mode_families
        pc%ranks_per_family(i) = mcpu
        if ((r.gt.0).and.(i.le.r))  pc%ranks_per_family(i) = pc%ranks_per_family(i) + 1 ! add extra rank if avaiable
      enddo
    else
      do i = 1, pc%n_mode_families
        pc%ranks_per_family(i) = ranks_per_family(i)
      enddo
    endif

    pc%rank_range(1) = 1
    do i = 2, pc%n_mode_families+1
      pc%rank_range(i) = pc%rank_range(i-1) + pc%ranks_per_family(i-1)
    enddo

    ! check for consistency
    r = 0
    do i= 2, pc%n_mode_families + 1
      r = r + pc%rank_range(i) - pc%rank_range(i-1)
    enddo
    if (r.ne.n_cpu) then
      write(*,*) "Error in distribution of ranks"
      call exit(0)
    endif

    do i = 1, n_cpu
      do j = 2, pc%n_mode_families + 1
        if ((i.ge.pc%rank_range(j-1)).and.(i.lt.pc%rank_range(j))) then
          rank_id(i) = j - 1
          exit
        endif
      enddo
    enddo

    do j=1,pc%n_mode_families
      r = 0
      do i = 1, n_cpu
        if (rank_id(i).eq.j) then
          r = r + 1
          pc%mode_families_ranks(j,r) = i - 1
        endif
      enddo
    enddo

    pc%family_id = rank_id(pc%my_id + 1)
    deallocate(rank_id)

    return

  end subroutine distribute_ranks
  
  !> Distribute RHS vector
  subroutine update_pc_rhs(pc,rhs_vec)
    use data_structure, only: type_PRECOND, type_RHS
    use mod_integer_types

    implicit none
    
    type(type_RHS)     :: rhs_vec
    type(type_PRECOND) :: pc

    integer(kind=int_all) :: i

    do i=1, pc%rhs%n
      pc%rhs%val(i) = rhs_vec%val(pc%row_index(i))
    enddo

  end subroutine update_pc_rhs
  
  !> Collect the solution vector from the individual mode groups
  subroutine gather_solution(pc,sol_vec)
    use data_structure, only: type_PRECOND, type_RHS
    use mod_integer_types
    use mpi_mod

    implicit none
    
    type(type_RHS)        :: sol_vec
    type(type_PRECOND)    :: pc

    integer               :: counts, ierr
    integer(kind=int_all) :: i  
    
    sol_vec%val(1:sol_vec%n) = 0.d0
    if (pc%my_id_n.eq.0) then
      do i = 1, pc%rhs%n
        sol_vec%val(pc%row_index(i)) = pc%rhs%val(i)*pc%row_factor
      enddo
    endif
    counts = sol_vec%n
    call MPI_AllReduce(MPI_IN_PLACE,sol_vec%val,counts,MPI_DOUBLE_PRECISION,MPI_SUM,pc%comm,ierr)
    
    return
  end subroutine gather_solution  

!> Deallocate arrays and reset to the default values
  subroutine reset_preconditioner(pc)
    use data_structure, only: type_PRECOND
    implicit none

    type(type_PRECOND) :: pc !, pc_def

    if (pc%initialized) then

      if (pc%structured) then

        deallocate(pc%rhs%val)
        pc%rhs%val => Null()
        deallocate(pc%row_index)
        pc%row_index => Null()
        deallocate(pc%send_counts, pc%recv_counts)
        pc%send_counts => Null()
        pc%recv_counts => Null()
        deallocate(pc%send_disp, pc%recv_disp)
        pc%send_disp => Null()
        pc%recv_disp => Null()
        deallocate(pc%istart, pc%ifinish)
        pc%istart => Null()
        pc%ifinish => Null()
        deallocate(pc%n_per_rank)
        pc%n_per_rank => Null()

        deallocate(pc%mat%val)
        deallocate(pc%mat%irn)
        deallocate(pc%mat%jcn)
        pc%mat%val => Null()
        pc%mat%irn => Null()
        pc%mat%jcn => Null()

        if (pc%mat%scaled) then
          deallocate(pc%mat%column_scaling)
          pc%mat%column_scaling => Null()
        endif

        pc%mat%scaled = .false.
        pc%mat%row_distributed = .false.
        pc%mat%col_distributed = .false.
        pc%mat%indexing = 1
        pc%mat%block_size = 1
        pc%structured = .false.

      endif

      deallocate(pc%mode_families_ranks)
      pc%mode_families_ranks => Null()
      deallocate(pc%mode_families_modes)
      pc%mode_families_modes => Null()
      deallocate(pc%mode_set)
      pc%mode_set => Null()

      pc%initialized = .false.

    endif

    return

  end subroutine reset_preconditioner


end module mod_preconditioner
