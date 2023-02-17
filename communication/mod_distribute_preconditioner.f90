module mod_distribute_preconditioner
#if !defined(DIRECT_CONSTRUCTION)
  use mod_integer_types

  implicit none

  private
  public update_pc_mat

contains

  !> Extract Preconditioner (PC) matrices from distributed global sparce matrix
  !! a_mat%val(1:nz_glob), a_mat%irn(1:nz_glob), a_mat%jcn(1:nz_glob)
  !! Uses splitted communication if number of send/recv entries exceeds INT_MAX
  !! nsplit - number of split communications
  !! nz_split - number of nonzeros to go through in each communication cycle
  !!
  !! Sends the reduced local matrices to the masters only
  !!  (centralize_harm_mat=.true.) or distribute by rows among all ranks
  !!
  !!    pc%mat%val(1:pc%mat%nnz), pc%rhs(1:pc%mat%n)
  !!    pc%mat%irn(1:pc%mat%nnz)
  !!    pc%mat%jcn(1:pc%mat%nnz)
  subroutine update_pc_mat(pc, a_mat, sim)

    use mod_parameters, only : n_tor, n_var
    use mpi_mod
    use mod_integer_types
    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS
    use mod_simulation_data, only: type_MHD_SIM    
    
    implicit none
    
    type(type_PRECOND)                 :: pc
    type(type_SP_MATRIX)               :: a_mat
    type(type_MHD_SIM), optional       :: sim

    integer                            :: my_id, n_cpu, j, k, l, n, ierr
    integer                            :: nm, ji, nr, lmode, kmode, n_i, n_j, isplit, icpu, ie
    integer(kind=int_all)              :: i, i0, i1, nz_split, ibufsize, block_size
    integer(kind=int_all)              :: n_tor_int

    real*8,  allocatable               :: Asnd_buffer(:)
    integer(kind=int_all), allocatable :: isnd_buffer(:), jsnd_buffer(:)
    real*8,  allocatable               :: Arcv_buffer(:)
    integer(kind=int_all), allocatable :: ircv_buffer(:), jrcv_buffer(:)
    integer(kind=int_all), allocatable :: indx(:)

    logical                            :: distribute_row, distribute_col

    integer :: cc, cr
    real t0, t1
    
    my_id   = pc%my_id
    n_cpu   = pc%n_cpu

    if (my_id .eq. 0) then
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*    distributing PC matrix     *'
      write(*,*) my_id,'*********************************'
    endif

    call system_clock(count=cc, count_rate=cr); t0 =  real(cc)/cr

! --- Copy of n_tor as long-integer for modulo functions (just to keep safe)
    n_tor_int = n_tor

    distribute_row = pc%mat%row_distributed.and.(.not.pc%mat%col_distributed)
    distribute_col = pc%mat%col_distributed.and.(.not.pc%mat%row_distributed)

    allocate(indx(n_cpu))

    if (.not.pc%structured) call set_pc_structure(pc,a_mat)

! --- Allocate PC matrices
    if (.not.associated(pc%mat%val)) allocate(pc%mat%val(pc%mat%nnz))
    if (.not.associated(pc%mat%irn)) allocate(pc%mat%irn(pc%mat%nnz))
    if (.not.associated(pc%mat%jcn)) allocate(pc%mat%jcn(pc%mat%nnz))

! --- Loop over communication splits
    do isplit = 1, pc%nsplit

      ibufsize = sum(pc%send_counts(isplit,1:n_cpu))

      allocate(Asnd_buffer(ibufsize))
      allocate(isnd_buffer(ibufsize))
      allocate(jsnd_buffer(ibufsize))

      call system_clock(count=cc, count_rate=cr); t1 =  real(cc)/cr

    ! prepare data to be distributed from the current rank
      indx(1) = 0 ! starting index for a particular destination rank
      do i = 2, n_cpu
        indx(i) = indx(i-1) + pc%send_counts(isplit,i-1)
      enddo

      if (pc%autodistribute_modes) then

        do i = pc%istart(isplit), pc%ifinish(isplit)
          n_i = (mod(a_mat%irn(i)-Int1,n_tor_int) + 1)/2
          n_j = (mod(a_mat%jcn(i)-Int1,n_tor_int) + 1)/2
          if (n_i .eq. n_j) then
            j = n_i + 1
            ji = pc%rank_range(j)
            if (distribute_row) then
              nr = pc%ranks_per_family(j)
              ji =  ji + min((a_mat%irn(i)-Int1)/pc%n_per_rank(j), nr-1) ! row bin index for j-th family
            elseif (distribute_col) then
              nr = pc%ranks_per_family(j)
              ji =  ji + min((a_mat%jcn(i)-Int1)/pc%n_per_rank(j), nr-1) ! column bin index for j-th family
            endif
            indx(ji) = indx(ji) + 1
            Asnd_buffer(indx(ji)) = a_mat%val(i)
            isnd_buffer(indx(ji)) = a_mat%irn(i)
            jsnd_buffer(indx(ji)) = a_mat%jcn(i)
          endif
        enddo

      else

  !$omp do private(i, j, ji, nm, nr, k, kmode, l, lmode, n_i, n_j)
        do i = pc%istart(isplit), pc%ifinish(isplit)
          n_i = mod(a_mat%irn(i)-Int1,n_tor_int) + 1
          n_j = mod(a_mat%jcn(i)-Int1,n_tor_int) + 1

          do j = 1, pc%n_mode_families
            nm = pc%modes_per_family(j) ! number of modes per j-th family
            do k = 1, nm
              kmode = pc%mode_families_modes(j,k)
              do l = 1, nm
                lmode = pc%mode_families_modes(j,l)
                if ((n_i.eq.kmode).and.(n_j.eq.lmode)) then
                  ji = pc%rank_range(j)
                  if (distribute_row) then
                    nr = pc%ranks_per_family(j)
                    ji =  ji + min((a_mat%irn(i)-1)/pc%n_per_rank(j), nr-1) ! row bin index for j-th family
                  elseif (distribute_col) then
                    nr = pc%ranks_per_family(j)
                    ji =  ji + min((a_mat%jcn(i)-1)/pc%n_per_rank(j), nr-1) ! row bin index for j-th family
                  endif
                  indx(ji) = indx(ji) + 1
                  Asnd_buffer(indx(ji)) = a_mat%val(i)
                  isnd_buffer(indx(ji)) = a_mat%irn(i)
                  jsnd_buffer(indx(ji)) = a_mat%jcn(i)
                endif
              enddo
            enddo
          enddo
        enddo

      endif

      allocate(Arcv_buffer(sum(pc%recv_counts(isplit,1:n_cpu))))
      call mpi_alltoallv(Asnd_buffer, pc%send_counts(isplit,1:n_cpu), pc%send_disp(isplit,1:n_cpu), MPI_DOUBLE_PRECISION, &
                         Arcv_buffer, pc%recv_counts(isplit,1:n_cpu), pc%splt_disp(isplit,1:n_cpu), MPI_DOUBLE_PRECISION,a_mat%comm,ierr)
      do icpu = 1, n_cpu
        i0 = pc%recv_disp(isplit,icpu) + int1
        i1 = i0 + pc%recv_counts(isplit,icpu) - int1
        pc%mat%val(i0:i1) = Arcv_buffer(pc%splt_disp(isplit,icpu) + 1:pc%splt_disp(isplit,icpu) + pc%recv_counts(isplit,icpu))
      enddo
      deallocate(Arcv_buffer)
      deallocate(Asnd_buffer)

      allocate(ircv_buffer(sum(pc%recv_counts(isplit,1:n_cpu))))
      call mpi_alltoallv(isnd_buffer, pc%send_counts(isplit,1:n_cpu), pc%send_disp(isplit,1:n_cpu), MPI_INTEGER_ALL, &
                         ircv_buffer, pc%recv_counts(isplit,1:n_cpu), pc%splt_disp(isplit,1:n_cpu), MPI_INTEGER_ALL,a_mat%comm,ierr)

      do icpu = 1, n_cpu
        i0 = pc%recv_disp(isplit,icpu) + int1
        i1 = i0 + pc%recv_counts(isplit,icpu) - int1
        pc%mat%irn(i0:i1) = ircv_buffer(pc%splt_disp(isplit,icpu) + 1:pc%splt_disp(isplit,icpu) + pc%recv_counts(isplit,icpu))
      enddo
      deallocate(ircv_buffer)
      deallocate(isnd_buffer)

      allocate(jrcv_buffer(sum(pc%recv_counts(isplit,1:n_cpu))))
      call mpi_alltoallv(jsnd_buffer, pc%send_counts(isplit,1:n_cpu), pc%send_disp(isplit,1:n_cpu), MPI_INTEGER_ALL, &
                         jrcv_buffer, pc%recv_counts(isplit,1:n_cpu), pc%splt_disp(isplit,1:n_cpu), MPI_INTEGER_ALL,a_mat%comm,ierr)
      do icpu = 1, n_cpu
        i0 = pc%recv_disp(isplit,icpu) + int1
        i1 = i0 + pc%recv_counts(isplit,icpu) - int1
        pc%mat%jcn(i0:i1) = jrcv_buffer(pc%splt_disp(isplit,icpu) + 1:pc%splt_disp(isplit,icpu) + pc%recv_counts(isplit,icpu))
      enddo
      deallocate(jrcv_buffer)
      deallocate(jsnd_buffer)

    enddo
    
    deallocate(indx)

! --- Change indices of the local matrices to local indices
!$omp do private(i,j,n_i,n_j)
    do i=1,pc%mat%nnz
      n_i = mod(pc%mat%irn(i)-Int1,n_tor_int) + 1
      do j=1, pc%mode_set_n
        if (n_i.eq.pc%mode_set(j)) then
          pc%mat%irn(i) = int((pc%mat%irn(i)-Int1)/n_tor_int)*pc%mode_set_n + j
          exit
        endif
      enddo

      n_j = mod(pc%mat%jcn(i)-Int1,n_tor_int) + 1
      do j=1, pc%mode_set_n
        if (n_j.eq.pc%mode_set(j)) then
          pc%mat%jcn(i) = int((pc%mat%jcn(i)-Int1)/n_tor_int)*pc%mode_set_n + j
          exit
        endif
      enddo
    enddo

    call system_clock(count=cc, count_rate=cr); t1 =  real(cc)/cr
!    write(*,*) "PC row range:", pc%my_id, minval(pc%mat%irn(1:pc%mat%nnz)), maxval(pc%mat%irn(1:pc%mat%nnz))
    if (pc%my_id.eq.0) write(*,'(A42,F6.2)') " Elapsed time updating preconditioner (s):",t1-t0
    return

  end subroutine update_pc_mat

  subroutine set_pc_structure(pc,a_mat)
    use mpi_mod
    use mod_integer_types
    use mod_parameters, only : n_tor, n_var
    use data_structure, only: type_SP_MATRIX, type_PRECOND    
    
    implicit none
    
    type(type_PRECOND)                 :: pc
    type(type_SP_MATRIX)               :: a_mat
    integer(kind=int_all), allocatable :: long_recv_counts(:), long_send_counts(:)
    integer(kind=int_all)              :: block_size, block_size2
    integer(kind=int_all)              :: nz_split, i0, i1, ind
    integer(kind=int_all)              :: nzg, nnz, nm
    integer                            :: isplit, i, j, ie, icpu
    integer                            :: ierr
    integer                            :: nr
    
    if (pc%my_id.eq.0) write(*,*) "Analyzing preconditioner"
    
    pc%mat%ng = (pc%mode_set_n)*(a_mat%ng)/n_tor ! rank of local PC matrix
    pc%mat%nr = pc%mat%ng
    pc%mat%nc = pc%mat%ng
    pc%n_glob = a_mat%ng
    
    allocate(pc%n_per_rank(pc%n_mode_families))
    
    ! if disributing, split the rows of the global matrix between ranks of families    
    do j = 1, pc%n_mode_families
      block_size = n_var*pc%modes_per_family(j)
      nr = pc%ranks_per_family(j)                                           ! number of ranks per j-th family
      pc%n_per_rank(j) = block_size*((a_mat%ng/block_size)/nr) - block_size ! number of rows per rank for j-th family (lower limit)
    enddo
 
    pc%mat%block_size = n_var*pc%mode_set_n ! set block size for current family

    allocate(long_send_counts(pc%n_cpu))
    allocate(long_recv_counts(pc%n_cpu))

    call MPI_Allreduce(a_mat%nnz,nzg,1,MPI_INTEGER_ALL,MPI_SUM,a_mat%comm,ierr)
    nm = sum(pc%modes_per_family(1:pc%n_mode_families)) ! number of modes in all families (can be larger than n_tor)
    nzg = (nzg/n_tor**2)*nm**2                          ! number of elements in all families (can be larger than global nnz)
    pc%nsplit = nzg/INT_MAX + 1

    if ((pc%my_id.eq.0).and.(pc%nsplit>1)) write(*,*) "Using split communication for preconditioner construction", pc%nsplit

    allocate(pc%istart(pc%nsplit),pc%ifinish(pc%nsplit))

    ! split global nz keeping it integer of (n_var*n_tor)**2
    block_size2 = (n_var*n_tor)**2
    nz_split = ((a_mat%nnz/block_size2)/pc%nsplit)*block_size2

    ! distribute indices for split communication
    pc%istart(1) = 1
    pc%ifinish(1) = nz_split
    
    do isplit = 2, pc%nsplit
      pc%istart(isplit)  = pc%ifinish(isplit-1) + 1
      pc%ifinish(isplit) = pc%istart(isplit) + nz_split - 1
    enddo
    pc%ifinish(pc%nsplit) = a_mat%nnz

    allocate(pc%send_counts(pc%nsplit,pc%n_cpu))
    allocate(pc%recv_counts(pc%nsplit,pc%n_cpu))
    allocate(pc%send_disp(pc%nsplit,pc%n_cpu))
    allocate(pc%recv_disp(pc%nsplit,pc%n_cpu))
    allocate(pc%splt_disp(pc%nsplit,pc%n_cpu))

    do isplit = 1, pc%nsplit
      i0 = pc%istart(isplit); i1 = pc%ifinish(isplit)
      call get_send_recv(i0,i1,long_send_counts,long_recv_counts,pc,a_mat)
      pc%send_counts(isplit,1:pc%n_cpu) = long_send_counts(1:pc%n_cpu)
      pc%recv_counts(isplit,1:pc%n_cpu) = long_recv_counts(1:pc%n_cpu)
    enddo

    nnz = 0
    do isplit = 1, pc%nsplit
      nnz = nnz + sum(pc%recv_counts(isplit,1:pc%n_cpu))
    enddo
    pc%mat%nnz=nnz

    do isplit = 1, pc%nsplit
      pc%send_disp(isplit,1) = 0
      do i = 2, pc%n_cpu
          pc%send_disp(isplit,i) = pc%send_disp(isplit,i-1) + pc%send_counts(isplit,i-1)
      enddo
    enddo

    !cpu:    |   cpu0    |   cpu1    |   cpu2    |
    !split:  | 1 | 2 | 3 | 1 | 2 | 3 | 1 | 2 | 3 |
    ! long displacements calculated cpu by cpu

    pc%recv_disp(1,1) = 0
    do isplit = 2, pc%nsplit
      pc%recv_disp(isplit,1) = pc%recv_disp(isplit-1,1) + pc%recv_counts(isplit-1,1)
    enddo
    do i = 2, pc%n_cpu
      pc%recv_disp(1,i) = pc%recv_disp(pc%nsplit,i-1) + pc%recv_counts(pc%nsplit,i-1)
      do isplit = 2, pc%nsplit
        pc%recv_disp(isplit,i) = pc%recv_disp(isplit-1,i) + pc%recv_counts(isplit-1,i)
      enddo
    enddo

    ! short int displacement for every isplit
    do isplit = 1, pc%nsplit
      pc%splt_disp(isplit,1) = 0
      do i = 2, pc%n_cpu
        pc%splt_disp(isplit,i) = pc%splt_disp(isplit,i-1) + pc%recv_counts(isplit,i-1)
      enddo
    enddo

    deallocate(long_send_counts,long_recv_counts)

    pc%mat%comm = pc%MPI_COMM_N    
    ! centralized
    allocate(pc%mat%index_min(pc%n_cpu_n)); pc%mat%index_min(1:pc%n_cpu_n) = 0
    allocate(pc%mat%index_max(pc%n_cpu_n)); pc%mat%index_max(1:pc%n_cpu_n) = 0
    pc%mat%index_min(pc%my_id_n + 1) = 1
    pc%mat%index_max(pc%my_id_n + 1) = pc%mat%nnz/(pc%mat%block_size*pc%mat%block_size)
    call MPI_Allreduce(MPI_IN_PLACE,pc%mat%index_min,pc%n_cpu_n,MPI_INTEGER_ALL,MPI_SUM,pc%MPI_COMM_N,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,pc%mat%index_max,pc%n_cpu_n,MPI_INTEGER_ALL,MPI_SUM,pc%MPI_COMM_N,ierr)
    
    allocate(pc%rhs%val(pc%mat%ng))
    pc%rhs%val(1:pc%mat%ng) = 0.d0
    pc%rhs%n = pc%mat%ng

    allocate(pc%row_index(pc%rhs%n))
    
    do ind = 0, a_mat%ng/n_tor - 1
      do i = 1, pc%mode_set_n
        pc%row_index(i + ind*pc%mode_set_n) =  pc%mode_set(i) + ind*n_tor
      enddo
    enddo
    
    write(*,*) pc%my_id, "PC matrix: ng, nnz", pc%mat%ng, pc%mat%nnz
    pc%structured = .true.

    return         
  
  end subroutine set_pc_structure

  !> Calculate send-recv counts
  subroutine get_send_recv(i0,i1,long_send_counts,long_recv_counts,pc,a_mat)

    use mpi_mod
    use mod_integer_types
    use mod_parameters, only : n_tor, n_var
    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS    
    
    implicit none
    
    type(type_PRECOND)                 :: pc
    type(type_SP_MATRIX)               :: a_mat    

    integer(kind=int_all), intent(in)                 :: i0, i1
    integer(kind=int_all), allocatable, intent(inout) :: long_recv_counts(:), long_send_counts(:)
    integer(kind=int_all), allocatable                :: sendrecv(:)
    integer(kind=int_all)                             :: i, n_tor_int
    integer                                           :: j, ji, nm, nr, k, kmode, l, lmode, n_i, n_j, ierr
    logical                                           :: distribute_row, distribute_col


    distribute_row = pc%mat%row_distributed.and.(.not.pc%mat%col_distributed)
    distribute_col = pc%mat%col_distributed.and.(.not.pc%mat%row_distributed)

    n_tor_int = n_tor

    long_send_counts(1:pc%n_cpu) = 0 ! number of elements to be sent from current rank to others

   ! calculate number of entries to be distributed from the current rank
    if (pc%autodistribute_modes) then

      do i = i0, i1
        n_i = (mod(a_mat%irn(i)-Int1,n_tor_int) + 1) / 2
        n_j = (mod(a_mat%jcn(i)-Int1,n_tor_int) + 1) / 2
        if (n_i .eq. n_j) then
          j = n_i + 1
          ji = pc%rank_range(j)
          if (distribute_row) then
            nr = pc%ranks_per_family(j)
            ji =  ji + min((a_mat%irn(i)-Int1)/pc%n_per_rank(j), nr-1) ! row bin index for j-th family
          elseif (distribute_col) then
            nr = pc%ranks_per_family(j)
            ji =  ji + min((a_mat%jcn(i)-Int1)/pc%n_per_rank(j), nr-1) ! column bin index for j-th family
          endif
          long_send_counts(ji) = long_send_counts(ji) + 1
        endif
      enddo

    else

!$omp do private(i, j, ji, nm, nr, k, kmode, l, lmode, n_i, n_j)
      do i=i0, i1
        n_i = mod(a_mat%irn(i)-Int1,n_tor_int) + 1
        n_j = mod(a_mat%jcn(i)-Int1,n_tor_int) + 1

        do j = 1, pc%n_mode_families
          nm = pc%modes_per_family(j) ! number of modes per j-th family
          do k = 1, nm
            kmode = pc%mode_families_modes(j,k)
            do l = 1, nm
              lmode = pc%mode_families_modes(j,l)
              if ((n_i.eq.kmode).and.(n_j.eq.lmode)) then
                ji = pc%rank_range(j)
                if (distribute_row) then
                  nr = pc%ranks_per_family(j)
                  ji =  ji + min((a_mat%irn(i)-Int1)/pc%n_per_rank(j), nr-1) ! row bin index for j-th family
                elseif (distribute_col) then
                  nr = pc%ranks_per_family(j)
                  ji =  ji + min((a_mat%jcn(i)-Int1)/pc%n_per_rank(j), nr-1) ! column bin index for j-th family
                endif
                long_send_counts(ji) = long_send_counts(ji) + 1
              endif
            enddo
          enddo
        enddo
      enddo

    endif

    allocate(sendrecv(pc%n_cpu*pc%n_cpu))
    sendrecv(1:pc%n_cpu*pc%n_cpu) = 0
    sendrecv(pc%my_id*pc%n_cpu + 1:(pc%my_id+1)*pc%n_cpu) = long_send_counts(1:pc%n_cpu)
    j = pc%n_cpu*pc%n_cpu
    call MPI_Barrier(a_mat%comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,sendrecv,j,MPI_INTEGER_ALL,MPI_SUM,a_mat%comm,ierr)

    long_recv_counts(1:pc%n_cpu) = 0
    do i = 1, pc%n_cpu
      long_recv_counts(i) = sendrecv(pc%n_cpu*(i-1) + pc%my_id + 1)
    enddo

    deallocate(sendrecv)

    return

  end subroutine get_send_recv
#endif
end module mod_distribute_preconditioner
