module strumpack_module
#ifdef USE_STRUMPACK

  use iso_c_binding
  use mpi
  use mod_integer_types

  implicit none
  type(C_PTR) :: spss ! STRUMPACK sparse solver
  integer(kind=C_INT_ALL), allocatable, target :: dist(:)  ! row distribution
  logical :: spss_initialized, spss_analyzed

  private
  public :: strumpack_init, strumpack_set_mat, strumpack_analyze, &
            strumpack_factorize, strumpack_solve, strumpack_finalize, &
            spss_initialized, spss_analyzed

  interface
    subroutine spk() bind(C)
      use iso_c_binding
    end subroutine spk

    subroutine spk_init(spss,comm) bind(C)
      use iso_c_binding
      implicit none

      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_init

    subroutine spk_set_mat(n,dist,irn,jcn,val,spss,comm,upd) bind(C)
      use iso_c_binding
      use mod_integer_types
      implicit none

      integer(kind=C_INT_ALL), dimension(:), pointer, intent(in) :: irn,jcn,dist
      real(kind=C_DOUBLE), dimension(:), pointer, intent(in) :: val
      integer(kind=C_INT_ALL), intent(in) :: n
      integer, intent(in) :: comm
      type(c_ptr), intent(inout) :: spss
      logical :: upd
    end subroutine spk_set_mat

    subroutine spk_reord(spss,comm) bind(C)
      use iso_c_binding
      implicit none

      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_reord

    subroutine spk_fact(spss,comm) bind(C)
      use iso_c_binding
      implicit none

      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_fact

    subroutine spk_solve(n,dist,rhsc,spss,comm) bind(C)
      use iso_c_binding
      use mod_integer_types
      implicit none

      integer(kind=C_INT_ALL), intent(in) :: n
      integer(kind=C_INT_ALL), dimension(:), pointer, intent(in) :: dist
      type(c_ptr), intent(inout) :: spss, rhsc
      integer, intent(in) :: comm
    end subroutine spk_solve

    subroutine spk_finalize(spss,comm) bind(C)
      use iso_c_binding
      implicit none

      type(c_ptr), intent(inout) :: spss
      integer, intent(in) :: comm
    end subroutine spk_finalize

  end interface

  contains
    subroutine strumpack_init(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr

        call spk_init(spss,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_init

    subroutine strumpack_set_mat(n,nnz,irn,jcn,val,block_size,comm,update,distributed,equilibrium) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        use sorting_module, only : remove_duplicates, convert2csr, convert_sorting
        use mod_integer_types

        implicit none

        integer comm,ierr
        integer(kind=C_INT_ALL), dimension(:), pointer :: irn, jcn, irnl, jcnl
        real(kind=C_DOUBLE),  dimension(:), pointer :: val, vall
        integer(kind=C_INT_ALL), intent(in) :: n
        integer(kind=C_INT_ALL), intent(inout) :: nnz
        integer, intent(in) :: block_size
        logical,intent(in),optional :: update, distributed, equilibrium

        integer(kind=C_INT_ALL), dimension(:), pointer :: myelm
        logical :: upd=.false., dflag=.false., eql=.false.

        integer :: rank, ncpu, indx=1
        integer(kind=int_all) :: nnzloc, nloc, i, j, imin, imax

        if(present(update)) upd = update
        if(present(distributed)) dflag = distributed
        if(present(equilibrium)) eql = equilibrium

        call MPI_COMM_RANK(comm, rank, ierr)
        call MPI_COMM_SIZE(comm, ncpu, ierr)

        if ((minval(irn).lt.1).or.(maxval(irn).gt.n)) then
          write(*,*) rank, ": Error: inconsistent row indices", minval(irn), maxval(irn), n
          call exit(0)
        endif

        if (eql) then
          call distribute_rows(n,1,block_size)
          dist(:) = dist(:) - indx

#if (defined(USEMKL))
          call convert2csr(indx,n,n,nnz,irn,jcn,val)
#else
          call remove_duplicates(n,nnz,irn,jcn,val)
          call convert2csr(indx,n,n,nnz,irn,jcn,val)
#endif
          call spk_set_mat(n,dist,irn,jcn,val,spss,comm,upd)

        else

          if ((.not. dflag).and.(ncpu.gt.1)) then
            ! distribute rows between ncpu
            call distribute_rows(n,ncpu,block_size)
            if (rank.eq.0) write(*,*) "Matrix is not row-distributed. Distributing now."

            allocate(myelm(nnz))
            j = 1
            do i=1, nnz
              if ((irn(i)>= dist(rank+1)).and.(irn(i)<= (dist(rank+2)-1))) then
                myelm(j) = i
                j = j + 1
              endif
            enddo

            nnzloc = j - 1
            nloc = dist(rank+2) - dist(rank+1)

            allocate(irnl(nnzloc), jcnl(nnzloc), vall(nnzloc))

            do i = 1, nnzloc
              irnl(i) = irn(myelm(i)) - dist(rank+1) + 1       ! irn starts from 1
              jcnl(i) = jcn(myelm(i))                          ! jcn remains the same
              vall(i) = val(myelm(i))
            enddo
            dist(:) = dist(:) - indx ! convert ot c-indexing
#if (defined(USEMKL))
            call convert2csr(indx,nloc,n,nnzloc,irnl,jcnl,vall)
#else
            call convert_sorting(nnzloc,irnl,jcnl,vall,block_size,indx)
#endif
            call spk_set_mat(nloc,dist,irnl,jcnl,vall,spss,comm,upd)
            deallocate(myelm,irnl,jcnl,vall)

          elseif (dflag.and.(ncpu.gt.1)) then
            ! get row distribution from irn in case of pre-distributed matrix
            if (allocated(dist)) deallocate(dist)
            allocate(dist(ncpu+1))
            dist(1:ncpu+1) = 0
            imin = minval(irn(1:nnz))
            imax = maxval(irn(1:nnz))
            dist(rank+1) = imin

            if (rank.eq.(ncpu-1)) dist(rank+2) = imax + 1
            call MPI_Allreduce(MPI_IN_PLACE,dist,ncpu+1,MPI_INTEGER,MPI_SUM,comm,ierr)

            ! check for consistency
            ierr = 0
            if ((dist(1).ne.1)) ierr = 1
            do i = 2, ncpu+1
              if (.not.(dist(i)>dist(i-1))) ierr = 1
            enddo

            if (ierr.ne.0) then
              write(*,*) "Error in harmonic matrix distribution"
              call exit(2)
            endif

            nloc = dist(rank+2) - dist(rank+1)
            
            irn(1:nnz) = irn(1:nnz) - imin + 1 ! irn starts from 1
            dist(1:ncpu+1) = dist(1:ncpu+1) - indx

#if (defined(USEMKL))
            call convert2csr(indx,nloc,n,nnz,irn,jcn,val)
#else
            call convert_sorting(nnz,irn,jcn,val,block_size,indx)
#endif
            call spk_set_mat(nloc,dist,irn,jcn,val,spss,comm,upd)

          else ! case of ncpu = 1
            call distribute_rows(n,1,1)
            dist(:) = dist(:) - indx

#if (defined(USEMKL))
            call convert2csr(indx,n,n,nnz,irn,jcn,val)
#else
            call convert_sorting(nnz,irn,jcn,val,block_size,indx)
#endif
            call spk_set_mat(n,dist,irn,jcn,val,spss,comm,upd)

          endif

        endif

        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_set_mat

    subroutine strumpack_analyze(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr

        call spk_reord(spss,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_analyze

    subroutine strumpack_factorize(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr

        call spk_fact(spss,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_factorize

    subroutine strumpack_solve(n,rhs,comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        use mod_integer_types
        implicit none

        real(kind=C_DOUBLE),  dimension(:), pointer :: rhs
        integer(kind=C_INT_ALL), intent(in) :: n
        integer, intent(in) :: comm

        integer :: rank, ncpu, ierr, nloc
        type(C_PTR) :: rhsc

        call MPI_COMM_RANK(comm, rank, ierr)
        call MPI_COMM_SIZE(comm, ncpu, ierr)

        rhsc = c_loc(rhs);

        call spk_solve(n,dist,rhsc,spss,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_solve

    subroutine strumpack_finalize(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr

        call spk_finalize(spss,comm)
        call MPI_Barrier(comm,ierr);

        return
    end subroutine strumpack_finalize

    subroutine distribute_rows(n,ncpu,block_size) bind(C)
      !> Distribute rows between members of MPI group

        use, intrinsic :: iso_c_binding
        use mpi
        use mod_integer_types
        implicit none

        integer, intent(in) :: ncpu, block_size
        integer(kind=int_all), intent(in) :: n
        integer(kind=C_INT_ALL), dimension(:), allocatable :: nr
        integer :: ierr, i

        if (allocated(dist)) deallocate(dist)
        allocate(dist(ncpu+1))
        allocate(nr(ncpu))

        nr = block_size*((n/block_size)/ncpu)
        nr(ncpu) = nr(ncpu) + (n - sum(nr))

        dist(1) = 1
        do i=1, ncpu
          dist(i+1)= dist(i) + nr(i)
        enddo

        deallocate(nr)

        return
    end subroutine distribute_rows

#endif
end module strumpack_module

