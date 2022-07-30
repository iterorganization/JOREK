module mod_distribute_preconditioner
  use mod_integer_types

  implicit none

  logical                            :: analyzed = .false.
  integer, allocatable :: send_counts(:,:), recv_counts(:,:), send_disp(:,:), recv_disp(:,:), rank_of(:)
  integer(kind=int_all), allocatable :: istart(:), ifinish(:)
  integer                            :: nsplit

  private
  public distribute_harmonics, distribute_vector

contains

  !> Extract Preconditioner (PC) matrices from distributed global sparce matrix
  !! A_glob(1:nz_glob), irn_glob(1:nz_glob), jcn_glob(1:nz_glob)
  !! Uses splitted communication if number of send/recv entries exceeds INT_MAX
  !! nsplit - number of split communications
  !! nz_split - number of nonzeros to go through in each communication cycle
  !!
  !! Sends the reduced local matrices to the masters only
  !!  (centralize_harm_mat=.true.) or distribute by rows among all ranks
  !!
  !!    mumps_par%A(1:mumps_par%nz), mumsp_par%rhs(1:mumps_par%n)
  !!    mumps_par%irn(1:mumps_par%nz)
  !!    mumps_par%jcn(1:mumps_par%nz)
  subroutine distribute_harmonics(my_id,my_id_n,n_cpu)

    use tr_module
    use mod_parameters, only : n_tor, n_var
    use global_distributed_matrix, only : irn_glob, jcn_glob, A_glob, rhs_glob, nz_glob, ndof_glob
    use mumps_module, only : mumps_par
    use mpi_mod
    use mod_integer_types
    use phys_module, only : centralize_harm_mat, modes_per_family, mode_families_modes, n_mode_families, &
                            ranks_per_family, autodistribute_modes
    use preconditioner_module, only: my_mode_set_n, my_mode_set, mode_families_ranks, rank_range

    implicit none

    integer                            :: my_id, my_id_n, n_cpu, j, k, l, n, ierr
    integer                            :: nm, ji, nr, lmode, kmode, n_i, n_j, isplit
    integer(kind=int_all)              :: i, i0, i1, n_tor_int, nz_split, ibufsize, block_size

    real*8,  allocatable               :: Asnd_buffer(:), Rsnd_buffer(:)
    integer(kind=int_all), allocatable :: isnd_buffer(:), jsnd_buffer(:)
    integer(kind=int_all), allocatable :: long_recv_counts(:), long_send_counts(:), indx(:), n_per_rank(:)

    logical                            :: distribute
    integer(kind=int_all)              :: INT_MAX
    integer(kind=int_all), parameter   :: Int1=1

    integer :: cc, cr
    real t0, t1

    if (my_id .eq. 0) then
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*      distributing matrix      *'
      write(*,*) my_id,'*********************************'
    endif

! --- integer limit for short integers (normally 2147483647)
#ifdef INTSIZE64
  ! --- Not sure why, but it seems MPI fails even with counters below the long-int limit
  ! --- Maybe MPI has some internal working arrays that need to be larger than the counters? half seems to work well...
  INT_MAX = 1000000000
#else
  ! --- If we're not using long-ints, then there is nothing to split anyway
  INT_MAX = 2147000000
#endif

    !call system_clock(count=cc, count_rate=cr); t0 =  real(cc)/cr

! --- Copy of n_tor as long-integer for modulo functions (just to keep safe)
    n_tor_int = n_tor
    distribute = .not.centralize_harm_mat

    allocate(indx(n_cpu))

    allocate(n_per_rank(n_mode_families))
    do j = 1, n_mode_families
      block_size = n_var*modes_per_family(j)
      nr = ranks_per_family(j) ! number of ranks per j-th family
      n_per_rank(j) = block_size*((ndof_glob/block_size)/nr) - block_size ! number of rows per rank for j-th family
    enddo

! --- Calculate send-recv counts for each communication split and store it for the future
    if (.not.analyzed) then

      allocate(rank_of(n_mode_families))
      rank_of = rank_range(1:n_mode_families) ! starting rank index for each family

      allocate(long_send_counts(n_cpu),long_recv_counts(n_cpu))

      i0 = Int1; i1 = nz_glob
      call get_send_recv(my_id,n_cpu,i0,i1,long_send_counts,long_recv_counts)

      mumps_par%nz = sum(long_recv_counts(1:n_cpu)) ! summing up all recieves

      nsplit = maxval((/maxval(long_recv_counts),maxval(long_send_counts)/))/INT_MAX + 1
      if ((my_id.eq.0).and.(nsplit>1)) write(*,*) "Using split communication for preconditioner construction", nsplit

      allocate(istart(nsplit),ifinish(nsplit))

      ! split global nz keeping it integer of (n_var*n_tor)**2
      block_size = (4*n_var*n_tor_int)**2
      nz_split = ((nz_glob/block_size)/nsplit)*block_size

      ! distribute indices for split communication
      istart(1) = 1
      ifinish(1) = nz_split
      do isplit = 2, nsplit
        istart(isplit) = ifinish(isplit-1) + 1
        ifinish(isplit) = istart(isplit) + nz_split - 1
      enddo
      ifinish(nsplit) = nz_glob

      allocate(send_counts(nsplit,n_cpu),recv_counts(nsplit,n_cpu))
      allocate(send_disp(nsplit,n_cpu),recv_disp(nsplit,n_cpu))

      if (nsplit.eq.1) then
        send_counts(1,1:n_cpu) = long_send_counts(1:n_cpu)
        recv_counts(1,1:n_cpu) = long_recv_counts(1:n_cpu)
      else
        do isplit = 1, nsplit
          i0 = istart(isplit); i1 = ifinish(isplit)
          call get_send_recv(my_id,n_cpu,i0,i1,long_send_counts,long_recv_counts)
          send_counts(isplit,1:n_cpu) = long_send_counts(1:n_cpu)
          recv_counts(isplit,1:n_cpu) = long_recv_counts(1:n_cpu)
        enddo
      endif

      do isplit = 1, nsplit
        send_disp(isplit,1) = 0
        do j = 2, n_cpu
            send_disp(isplit,j) = send_disp(isplit,j-1) + send_counts(isplit,j-1)
        enddo
      enddo

      recv_disp(1,1) = 0
      do j = 2, n_cpu
          recv_disp(1,j) = recv_disp(1,j-1) + recv_counts(1,j-1)
      enddo
      do isplit = 2, nsplit
        recv_disp(isplit,1) = recv_disp(isplit-1,n_cpu) + recv_counts(isplit-1,n_cpu)
        do j = 2, n_cpu
            recv_disp(isplit,j) = recv_disp(isplit,j-1) + recv_counts(isplit,j-1)
        enddo
      enddo

      analyzed = .true.
      deallocate(long_send_counts,long_recv_counts)

    endif

! --- Allocate PC matrices
    if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"dh_mumps_par%A",CAT_DMATRIX)
    if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"dh_mumps_par%irn",CAT_DMATRIX)
    if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"dh_mumps_par%jcn",CAT_DMATRIX)

    call tr_allocatep(mumps_par%A,  Int1,mumps_par%nz,"dh_mumps_par%A",CAT_DMATRIX)
    call tr_allocatep(mumps_par%irn,Int1,mumps_par%nz,"dh_mumps_par%irn",CAT_DMATRIX)
    call tr_allocatep(mumps_par%jcn,Int1,mumps_par%nz,"dh_mumps_par%jcn",CAT_DMATRIX)

    mumps_par%A = 0.d0
    mumps_par%irn = 0
    mumps_par%jcn = 0

! --- Loop over communication splits
    do isplit = 1, nsplit

      ibufsize = sum(send_counts(isplit,1:n_cpu))

      call tr_allocate(Asnd_buffer,Int1,ibufsize,"dh_Asnd_buffer",CAT_DMATRIX)
      call tr_allocate(isnd_buffer,Int1,ibufsize,"dh_isnd_buffer",CAT_DMATRIX)
      call tr_allocate(jsnd_buffer,Int1,ibufsize,"dh_jsnd_buffer",CAT_DMATRIX)

      call system_clock(count=cc, count_rate=cr); t1 =  real(cc)/cr

    ! prepare data to be distributed from the current rank
      indx(1) = 0 ! starting index for a particular destination rank
      do i = 2, n_cpu
        indx(i) = indx(i-1) + send_counts(isplit,i-1)
      enddo

      if (autodistribute_modes) then

        do i = istart(isplit), ifinish(isplit)
          n_i = (mod(irn_glob(i)-Int1,n_tor_int) + 1) / 2
          n_j = (mod(jcn_glob(i)-Int1,n_tor_int) + 1) / 2
          if (n_i .eq. n_j) then
            j = n_i + 1
            ji = rank_of(j)
            if (distribute) then
              nr = ranks_per_family(j)
#ifdef USE_STRUMPACK
              ji =  ji + min((irn_glob(i)-Int1)/n_per_rank(j), nr-1) ! row bin index for j-th family
#elif USE_PASTIX6
              ji =  ji + min((jcn_glob(i)-Int1)/n_per_rank(j), nr-1) ! column bin index for j-th family
#endif
            endif
            indx(ji) = indx(ji) + 1
            Asnd_buffer(indx(ji)) = A_glob(i)
            isnd_buffer(indx(ji)) = irn_glob(i)
            jsnd_buffer(indx(ji)) = jcn_glob(i)
          endif
        enddo

      else

  !$omp do private(i, j, ji, nm, nr, k, kmode, l, lmode, n_i, n_j)
        do i = istart(isplit), ifinish(isplit)
          n_i = mod(irn_glob(i)-Int1,n_tor_int) + 1
          n_j = mod(jcn_glob(i)-Int1,n_tor_int) + 1

          do j = 1, n_mode_families
            nm = modes_per_family(j) ! number of modes per j-th family
            do k = 1, nm
              kmode = mode_families_modes(j,k)
              do l = 1, nm
                lmode = mode_families_modes(j,l)
                if ((n_i.eq.kmode).and.(n_j.eq.lmode)) then
                  ji = rank_of(j)
                  if (distribute) then
                    nr = ranks_per_family(j)
#ifdef USE_STRUMPACK
                    ji =  ji + min((irn_glob(i)-1)/n_per_rank(j), nr-1) ! row bin index for j-th family
#elif USE_PASTIX6
                    ji =  ji + min((jcn_glob(i)-1)/n_per_rank(j), nr-1) ! row bin index for j-th family
#endif
                  endif
                  indx(ji) = indx(ji) + 1
                  Asnd_buffer(indx(ji)) = A_glob(i)
                  isnd_buffer(indx(ji)) = irn_glob(i)
                  jsnd_buffer(indx(ji)) = jcn_glob(i)
                endif
              enddo
            enddo
          enddo
        enddo

      endif

      call mpi_alltoallv(Asnd_buffer,send_counts(isplit,1:n_cpu),send_disp(isplit,1:n_cpu),MPI_DOUBLE_PRECISION, &
                       mumps_par%A,  recv_counts(isplit,1:n_cpu),recv_disp(isplit,1:n_cpu),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call mpi_alltoallv(isnd_buffer,send_counts(isplit,1:n_cpu),send_disp(isplit,1:n_cpu),MPI_INTEGER_ALL, &
                       mumps_par%irn,recv_counts(isplit,1:n_cpu),recv_disp(isplit,1:n_cpu),MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)
      call mpi_alltoallv(jsnd_buffer,send_counts(isplit,1:n_cpu),send_disp(isplit,1:n_cpu),MPI_INTEGER_ALL, &
                       mumps_par%jcn,recv_counts(isplit,1:n_cpu),recv_disp(isplit,1:n_cpu),MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

      call tr_deallocate(Asnd_buffer,"dh_Asnd_buffer",CAT_DMATRIX)
      call tr_deallocate(isnd_buffer,"dh_isnd_buffer",CAT_DMATRIX)
      call tr_deallocate(jsnd_buffer,"dh_jsnd_buffer",CAT_DMATRIX)

    enddo

! --- Change indices of the local matrices to local indices
!$omp do private(i,j,n_i,n_j)
    do i=1,mumps_par%nz
      n_i = mod(mumps_par%irn(i)-Int1,n_tor_int) + 1
      do j=1, my_mode_set_n
        if (n_i.eq.my_mode_set(j)) then
          mumps_par%irn(i) = int((mumps_par%irn(i)-Int1)/n_tor_int)*my_mode_set_n + j
          exit
        endif
      enddo

      n_j = mod(mumps_par%jcn(i)-Int1,n_tor_int) + 1
      do j=1, my_mode_set_n
        if (n_j.eq.my_mode_set(j)) then
          mumps_par%jcn(i) = int((mumps_par%jcn(i)-Int1)/n_tor_int)*my_mode_set_n + j
          exit
        endif
      enddo
    enddo

    !call system_clock(count=cc, count_rate=cr); t1 =  real(cc)/cr
    !if (my_id.eq.0) write(*,*) "Elapsed time distributing (total):",t1-t0


    mumps_par%n =  my_mode_set_n*ndof_glob/n_tor

    if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"dh_mumps_par%rhs",CAT_DMATRIX)
    call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n,"dh_mumps_par%rhs",CAT_DMATRIX)

  end subroutine distribute_harmonics


  !> Distribute vector rhs among each mode group master
  subroutine distribute_vector(rhs,rhs_dis,comm)

    use global_distributed_matrix, only: ndof_glob
    use mumps_module, only: mumps_par
    use mpi_mod
    use preconditioner_module, only: my_row_index

    implicit none

    real*8                :: rhs(:), rhs_dis(:)
    integer               :: comm,ierr
    integer(kind=int_all) :: i

    call MPI_BCAST(rhs,ndof_glob,MPI_DOUBLE_PRECISION,0,comm,ierr)

    do i=1, mumps_par%n
      rhs_dis(i) = rhs(my_row_index(i))
    enddo

  end subroutine distribute_vector


  !> Calculate send-recv counts
  subroutine get_send_recv(my_id,n_cpu,i0,i1,long_send_counts,long_recv_counts)

    use mpi_mod
    use mod_integer_types
    use mod_parameters, only : n_tor, n_var
    use global_distributed_matrix, only : irn_glob, jcn_glob, A_glob, ndof_glob
    use phys_module, only : centralize_harm_mat, modes_per_family, mode_families_modes, n_mode_families, &
                            ranks_per_family, autodistribute_modes
    use preconditioner_module, only: mode_families_ranks, rank_range

    implicit none

    integer, intent(in) :: my_id, n_cpu
    integer(kind=int_all), intent(in)                 :: i0, i1
    integer(kind=int_all), allocatable, intent(inout) :: long_recv_counts(:), long_send_counts(:)
    integer(kind=int_all), allocatable                :: sendrecv(:), n_per_rank(:)
    integer(kind=int_all)                             :: i, n_tor_int, block_size
    integer(kind=int_all), parameter                  :: Int1=1
    integer                                           :: j, ji, nm, nr, k, kmode, l, lmode, n_i, n_j, ierr
    logical                                           :: distribute

    distribute = .not.centralize_harm_mat
    n_tor_int = n_tor
    allocate(n_per_rank(n_mode_families))
    do j = 1, n_mode_families
      block_size = n_var*modes_per_family(j)
      nr = ranks_per_family(j) ! number of ranks per j-th family
      n_per_rank(j) = block_size*((ndof_glob/block_size)/nr) - block_size ! number of rows per rank for j-th family
    enddo

    long_send_counts = 0 ! number of elements to be sent from current rank to others

   ! calculate number of entries to be distributed from the current rank
    if (autodistribute_modes) then

      do i = i0, i1
        n_i = (mod(irn_glob(i)-Int1,n_tor_int) + 1) / 2
        n_j = (mod(jcn_glob(i)-Int1,n_tor_int) + 1) / 2
        if (n_i .eq. n_j) then
          j = n_i + 1
          ji = rank_of(j)
          if (distribute) then
            nr = ranks_per_family(j)
#ifdef USE_STRUMPACK
            ji =  ji + min((irn_glob(i)-Int1)/n_per_rank(j), nr-1) ! row bin index for j-th family
#elif USE_PASTIX6
            ji =  ji + min((jcn_glob(i)-Int1)/n_per_rank(j), nr-1) ! column bin index for j-th family
#endif
          endif
          long_send_counts(ji) = long_send_counts(ji) + 1
        endif
      enddo

    else

!$omp do private(i, j, ji, nm, nr, k, kmode, l, lmode, n_i, n_j)
      do i=i0, i1
        n_i = mod(irn_glob(i)-Int1,n_tor_int) + 1
        n_j = mod(jcn_glob(i)-Int1,n_tor_int) + 1

        do j = 1, n_mode_families
          nm = modes_per_family(j) ! number of modes per j-th family
          do k = 1, nm
            kmode = mode_families_modes(j,k)
            do l = 1, nm
              lmode = mode_families_modes(j,l)
              if ((n_i.eq.kmode).and.(n_j.eq.lmode)) then
                ji = rank_of(j)
                if (distribute) then
                  nr = ranks_per_family(j)
#ifdef USE_STRUMPACK
                  ji =  ji + min((irn_glob(i)-Int1)/n_per_rank(j), nr-1) ! row bin index for j-th family
#elif USE_PASTIX6
                  ji =  ji + min((jcn_glob(i)-Int1)/n_per_rank(j), nr-1) ! column bin index for j-th family
#endif
                endif
                long_send_counts(ji) = long_send_counts(ji) + 1
              endif
            enddo
          enddo
        enddo
      enddo

    endif

    allocate(sendrecv(n_cpu*n_cpu))
    sendrecv = 0
    sendrecv(my_id*n_cpu + 1:(my_id+1)*n_cpu) = long_send_counts(1:n_cpu)
    j = n_cpu*n_cpu
    call mpi_allreduce(MPI_IN_PLACE,sendrecv,j,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)

    long_recv_counts = 0
    do i = 1, n_cpu
      long_recv_counts(i) = sendrecv(n_cpu*(i-1) + my_id + 1)
    enddo

    deallocate(sendrecv)

    return

  end subroutine get_send_recv

end module mod_distribute_preconditioner
