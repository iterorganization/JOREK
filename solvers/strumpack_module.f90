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

      integer(kind=C_INT_ALL), intent(in) :: n
      integer, intent(in) :: comm
      type(c_ptr) :: irn, jcn, val, dist
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

    subroutine spk_solve(n,dist,rhs,spss,comm) bind(C)
      use iso_c_binding
      use mod_integer_types
      implicit none

      integer(kind=C_INT_ALL), intent(in) :: n
      integer(kind=C_INT_ALL), dimension(:), pointer, intent(in) :: dist
      type(c_ptr), intent(inout) :: spss, rhs
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
    subroutine strumpack_init(comm)
        implicit none

        integer comm,ierr

        call spk_init(spss,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_init

    subroutine strumpack_set_mat(n,nnz,irn,jcn,val,block_size,comm,update,distributed,equilibrium)
        use, intrinsic :: iso_c_binding
        use sorting_module, only : remove_duplicates, convert2csr, convert_sorting
        use mod_integer_types

        implicit none

        integer comm,ierr
        integer(kind=C_INT_ALL), dimension(:), pointer :: irn, jcn, irn_d, jcn_d
        real(kind=C_DOUBLE),  dimension(:), pointer :: val, val_d

        integer(kind=C_INT_ALL), intent(in) :: n
        integer(kind=C_INT_ALL), intent(inout) :: nnz
        integer, intent(in) :: block_size
        logical,intent(in),optional :: update, distributed, equilibrium
        
        integer :: rank, ncpu
        integer(kind=int_all) :: nnz_d, n_d, i, j, imin, imax, indx        

        integer(kind=C_INT_ALL), dimension(:), pointer :: myelm
        logical :: upd, dflag, eql
        type(c_ptr) :: irn_c, jcn_c, val_c, dist_c
        
        upd = .false.
        dflag = .false.
        eql = .false.

        if(present(update)) upd = update
        if(present(distributed)) dflag = distributed
        if(present(equilibrium)) eql = equilibrium

        call MPI_COMM_RANK(comm, rank, ierr)
        call MPI_COMM_SIZE(comm, ncpu, ierr)

        indx = 1

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

          nnz_d = j - 1
          n_d = dist(rank+2) - dist(rank+1) ! number of local rows

          allocate(irn_d(nnz_d), jcn_d(nnz_d), val_d(nnz_d))

          do i = 1, nnz_d
            irn_d(i) = irn(myelm(i)) - dist(rank+1) + indx    ! irn starts from index
            jcn_d(i) = jcn(myelm(i))                          ! jcn remains the same
            val_d(i) = val(myelm(i))
          enddo

          deallocate(irn,jcn,val)
          irn => irn_d
          jcn => jcn_d
          val => val_d

          deallocate(myelm)

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

          n_d = dist(rank+2) - dist(rank+1)
          nnz_d = nnz
          irn(1:nnz_d) = irn(1:nnz_d) - imin + indx ! irn starts from indx
          
        elseif (ncpu.eq.1) then        
          call distribute_rows(n,1,block_size)
          n_d = n
          nnz_d = nnz    

        endif
        
        if (eql) then
          call remove_duplicates(n,nnz,irn,jcn,val)
          n_d = n
          nnz_d = nnz
        endif
        
        if (indx.eq.1) then
          irn(1:nnz_d) = irn(1:nnz_d) - indx;
          jcn(1:nnz_d) = jcn(1:nnz_d) - indx;
          dist(1:ncpu+1) = dist(1:ncpu+1) - indx
          indx = 0
        endif
        
#if (defined(USEMKL))
        irn_c = c_loc(irn); jcn_c = c_loc(jcn); val_c = c_loc(val); dist_c = c_loc(dist)
        call convert2csr(indx,n_d,n,nnz_d,irn_c,jcn_c,val_c)
#else
        call convert_sorting(nnz_d,irn,jcn,val,block_size,indx)
        irn_c = c_loc(irn); jcn_c = c_loc(jcn); val_c = c_loc(val); dist_c = c_loc(dist)
#endif

        call spk_set_mat(n_d,dist_c,irn_c,jcn_c,val_c,spss,comm,upd)

        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_set_mat

    subroutine strumpack_analyze(comm)
        implicit none

        integer comm,ierr

        call spk_reord(spss,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_analyze

    subroutine strumpack_factorize(comm) bind(C)
        implicit none

        integer comm,ierr

        call spk_fact(spss,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_factorize

    subroutine strumpack_solve(n,rhs,comm)
        use, intrinsic :: iso_c_binding
        use mod_integer_types
        implicit none

        real(kind=C_DOUBLE),  dimension(:), pointer :: rhs
        integer(kind=C_INT_ALL), intent(in) :: n
        integer, intent(in) :: comm

        integer :: rank, ncpu, ierr, n_d
        type(C_PTR) :: rhs_c

        call MPI_COMM_RANK(comm, rank, ierr)
        call MPI_COMM_SIZE(comm, ncpu, ierr)

        rhs_c = c_loc(rhs);

        call spk_solve(n,dist,rhs_c,spss,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine strumpack_solve

    subroutine strumpack_finalize(comm)
        implicit none

        integer comm,ierr
        
        if (allocated(dist)) deallocate(dist)

        call spk_finalize(spss,comm)
        call MPI_Barrier(comm,ierr);

        return
    end subroutine strumpack_finalize

    subroutine distribute_rows(n,ncpu,block_size)
      !> Distribute rows between members of MPI group
        use, intrinsic :: iso_c_binding
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

