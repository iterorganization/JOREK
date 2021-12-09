module mod_pastix
#ifdef USE_PASTIX6

  use iso_c_binding
  use mpi

  implicit none

  type(C_PTR) :: spm, pastix_data ! sparse solver (distributed)
  type(C_PTR) :: iparm, dparm

  integer(kind=C_INT), dimension(:), pointer :: loc2glob, glob2loc  ! mapping for column distribution
  real(kind=C_DOUBLE), allocatable   :: col_scaling(:)
  logical :: spm_initialized, spm_analyzed, spm_scaled

  private
  public :: pastix_init, pastix_set_mat, pastix_analyze, &
            pastix_factorize, pastix_solve, get_residual, &
            pastix_finalize, spm_initialized, spm_analyzed

  interface
    subroutine ptx() bind(C)
      use iso_c_binding
    end subroutine ptx

    subroutine ptx_init(pastix_data,spm,iparm,dparm,comm) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(out) :: pastix_data, spm
      type(C_PTR), intent(out) :: iparm, dparm
      integer, intent(in) :: comm

    end subroutine ptx_init

    subroutine ptx_set_mat(spm, indx, n, nnz, n_d, nnz_d, dof, rptr, cptr, values, &
               loc2glob, glob2loc, comm, update, check) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm
      integer(kind=C_INT) :: indx, comm
      integer(kind=C_INT) :: n, n_d, nnz, nnz_d, dof

      !type(C_PTR), intent(inout) :: rptr, cptr, values, loc2glob

      integer(kind=C_INT), dimension(:), pointer :: rptr, cptr, loc2glob, glob2loc
      real(kind=C_DOUBLE), dimension(:), pointer :: values
      logical :: update, check

    end subroutine ptx_set_mat

    subroutine ptx_analyze(pastix_data, spm) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm, pastix_data

    end subroutine ptx_analyze

    subroutine ptx_factorize(pastix_data, spm) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm, pastix_data

    end subroutine ptx_factorize

    subroutine ptx_solve(pastix_data, spm, rhsc, refine) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm, pastix_data
      type(C_PTR), intent(inout) :: rhsc
      logical :: refine

    end subroutine ptx_solve

    subroutine ptx_finalize(pastix_data, spm, iparm, dparm) bind(C)

      use iso_c_binding
      implicit none

      type(C_PTR), intent(inout) :: spm, pastix_data
      type(C_PTR), intent(inout) :: iparm, dparm

    end subroutine ptx_finalize

    subroutine get_residual(n, nnz, irn, jcn, val, x, rhs, indexing) bind(C)

      use iso_c_binding
      implicit none

      integer(kind=C_INT), dimension(:), pointer, intent(in) :: irn, jcn
      real(kind=C_DOUBLE), dimension(:), pointer, intent(in) :: val, rhs, x
      integer(kind=C_INT) :: n, nnz, indexing

    end subroutine get_residual

  end interface

  contains
!> Initialize PaStiX solver instance
    subroutine pastix_init(comm) bind(C)
        use, intrinsic :: iso_c_binding
        use mpi
        implicit none

        integer comm,ierr

        call ptx_init(pastix_data,spm,iparm,dparm,comm)
        call MPI_Barrier(comm,ierr)

        return
    end subroutine pastix_init

!> Prepare sparse matrix for pastix solver
!! Matrix is distributed column-wise among MPI processes in comm
!! The values are scaled such that the largest value in each column is 1
!! The matrix is converted to CSC format as required by distributed PaStiX
    subroutine pastix_set_mat(n, nnz, irn, jcn, val, block_size, comm, &
                update,distributed,equilibrium) bind(C)

      use, intrinsic :: iso_c_binding
      use mod_coicsr, only: coicsr, coicsr2
      use sorting_module, only: remove_duplicates, convert2csr
      implicit none

      integer(kind=C_INT), dimension(:), pointer :: irn, jcn
      real(kind=C_DOUBLE), dimension(:), pointer :: val

      logical,intent(in),optional :: update, distributed, equilibrium

      integer(kind=C_INT) :: indexing=1, block_size, block_size2, dof
      integer(kind=C_INT) :: n, n_d, nnz, nnz_d, jmin, jmax, i, n_block, nnz_block

      integer :: my_id_n, n_cpu_n, ierr, comm

      integer(kind=C_INT), dimension(:), allocatable :: iwk

      logical :: upd, dflag, eql
      
      upd = .false.
      dflag = .false.
      eql = .false.
      spm_scaled = .false.

      if(present(update)) upd = update
      if(present(distributed)) dflag = distributed
      if(present(equilibrium)) eql = equilibrium
      
      call MPI_Comm_rank(comm, my_id_n, ierr)

      ! if not already distributed distribute matrix column-wise
      if (.not.dflag) then
        call distribute_matrix(indexing, block_size, n, nnz, jmin, jmax, nnz_d, jcn, irn, val, comm)
      elseif (dflag) then
      ! already distributed;
        jmin = minval(jcn(1:nnz))
        jmax = maxval(jcn(1:nnz))
        nnz_d = nnz
        jcn(1:nnz_d) = jcn(1:nnz_d) - jmin + indexing
      endif
      
      if (.not.eql) then
        call do_column_scaling(n,nnz_d,irn,jcn,val,comm)
        spm_scaled = .true.
      endif      
      
      dof = 1
!#if (defined(USE_BLOCK))
!      if (.not.eql) then
!        n_block   = n/block_size
!        block_size2 = block_size*block_size
!        nnz_block = nnz_d/block_size2
!        if (block_size>1) then
!          do i=1,nnz_block
!            irn(i) = (irn((i-1)*block_size2 + 1) - 1)/block_size + 1
!            jcn(i) = (jcn((i-1)*block_size2 + 1) - 1)/block_size + 1
!          enddo
!        endif
!        jmin = minval(jcn(1:nnz_block))
!        jmax = maxval(jcn(1:nnz_block))
!        n = n_block
!        dof = block_size
!      endif
!#endif      
      
      n_d = jmax - jmin + 1
      allocate(loc2glob(n_d))
      do i = 1, n_d
        loc2glob(i) = i - 1 + jmin;
      enddo

      allocate(glob2loc(n)); glob2loc(1:n) = 0

      do i = 1, n_d
        glob2loc(loc2glob(i) + 1 - indexing)= - my_id_n - 1;
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, glob2loc, n, MPI_INTEGER, MPI_SUM, comm, ierr)

      do i = 1, n_d
        glob2loc(loc2glob(i) + 1 - indexing)= loc2glob(i) - loc2glob(1)
      enddo    

      if (eql) then
#if (defined(USEMKL))
        call convert2csr(indexing, n_d, n, nnz_d, jcn, irn, val)
#else
        call remove_duplicates(n, nnz_d, jcn, irn, val)
        call convert2csr(indexing, n_d, n, nnz_d, jcn, irn, val)
#endif
      else
        if (allocated(iwk)) deallocate(iwk)
        allocate(iwk(n + 1))
!#if (defined(USE_BLOCK))
!        call coicsr2(n,nnz_block,val,irn(1:nnz_block),jcn(1:nnz_block),block_size,iwk)
!#else
        call coicsr(n,nnz_d,1,val,irn,jcn,iwk)
!#endif
        deallocate(iwk)
      endif

      call MPI_Barrier(comm,ierr)

      call ptx_set_mat(spm, indexing, n, nnz, n_d, nnz_d, dof, irn, jcn, val, loc2glob, glob2loc, comm, upd, eql)

      return

    end subroutine pastix_set_mat

!> Perform matrix analysis/reordering
    subroutine pastix_analyze() bind(C)

      use iso_c_binding
      implicit none

      call ptx_analyze(pastix_data, spm)

    end subroutine pastix_analyze

!> Perform numerical LU factorization
    subroutine pastix_factorize() bind(C)

      use iso_c_binding
      implicit none

      call ptx_factorize(pastix_data, spm)

    end subroutine pastix_factorize

!> Calculate the solution
!! solution is placed into the rhs
    subroutine pastix_solve(n, rhs, refine) bind(C)

      use iso_c_binding
      implicit none

      real(kind=C_DOUBLE), pointer, intent(inout) :: rhs(:)
      logical, intent(in), optional :: refine
      type(C_PTR) :: rhsc
      integer, intent(in) :: n
      integer :: i

      logical :: ref
      ref = .false.

      if(present(refine)) ref = refine

      rhsc = c_loc(rhs)

      call ptx_solve(pastix_data, spm, rhsc, ref)
      
      if (spm_scaled) then
        do i = 1, n
          rhs(i) =  rhs(i)/col_scaling(i)
        enddo
      endif

    end subroutine pastix_solve

!> Finalize PaStiX solver instance
    subroutine pastix_finalize() bind(C)
      use iso_c_binding
      implicit none

      call ptx_finalize(pastix_data, spm,  iparm, dparm)
      deallocate(loc2glob, glob2loc)
      if (allocated(col_scaling))  deallocate(col_scaling)

    end subroutine pastix_finalize

!> Distribute matrix by the first array (irn or jcn) among MPI processes in comm
    subroutine distribute_matrix(indexing,block_size,n,nnz,imin,imax,nnz_d,irn,jcn,val,comm)
      use iso_c_binding
      implicit none

      integer(kind=C_INT), dimension(:), pointer :: irn, jcn
      real(kind=C_DOUBLE), dimension(:), pointer :: val

      integer, intent(in):: indexing, block_size
      integer, intent(in) :: n, nnz
      integer, intent(out) :: nnz_d, imin, imax

      integer(kind=C_INT), dimension(:), pointer :: irn_d, jcn_d
      real(kind=C_DOUBLE), dimension(:), pointer :: val_d
      integer :: n_cpu_n, my_id_n, comm, ierr
      integer :: i, j, indx

      integer(kind=C_INT), allocatable :: dist(:), myelm(:)

      call MPI_Comm_rank(comm, my_id_n, ierr)
      call MPI_Comm_size(comm, n_cpu_n, ierr)

      if (allocated(dist)) deallocate(dist)
      allocate(dist(n_cpu_n+1))

      ! number of rows/columns per cpu with the last one getting extra
      dist(2:n_cpu_n+1) = block_size*((n/block_size)/n_cpu_n)
      dist(n_cpu_n+1) = dist(n_cpu_n+1) + (n - sum(dist(2:n_cpu_n+1)))

      dist(1) = indexing
      do i=2, n_cpu_n+1
        dist(i) = dist(i) + dist(i-1)
      enddo
      
      imin = dist(my_id_n + 1)
      imax = dist(my_id_n + 2) - 1

      allocate(myelm(nnz))
      j = 1
      do i=1, nnz
        if ((irn(i)>= imin).and.(irn(i)<=imax)) then
          myelm(j) = i
          j = j + 1
        endif
      enddo

      nnz_d = j - 1

      allocate(irn_d(nnz_d),jcn_d(nnz_d),val_d(nnz_d))

      do i = 1, nnz_d
        irn_d(i) = irn(myelm(i)) - dist(my_id_n + 1) + indexing
        jcn_d(i) = jcn(myelm(i))
        val_d(i) = val(myelm(i))
      enddo

      deallocate(myelm,dist)
      deallocate(irn,jcn,val)

      irn => irn_d
      jcn => jcn_d
      val => val_d
      return

    end subroutine distribute_matrix
    
!> Perform column scaling after column-wise distribution
!! n - global number of rows
!! nnz - number of nnz in the local set of columns
!! the column_scaling array is defined externally
    subroutine do_column_scaling(n,nnz,irn,jcn,val,comm)
      use mod_integer_types
      use iso_c_binding
      use tr_module
  
      implicit none
      
      integer(kind=int_all) :: n, nnz
      integer(kind=C_INT), dimension(:), pointer :: irn, jcn
      real(kind=C_DOUBLE), dimension(:), pointer :: val
      
      integer(kind=int_all) :: i, j, k, jmin, jmax, n_d
      integer(kind=int_all), allocatable :: rcounts(:), displs(:)
      real(kind=C_DOUBLE), allocatable   :: loc_col_scaling(:)
      integer :: n_cpu_n, my_id_n, comm, ierr
      integer(kind=int_all), parameter   :: Int1=1
      
      call MPI_Comm_rank(comm, my_id_n, ierr)
      call MPI_Comm_size(comm, n_cpu_n, ierr)
    
    ! calculate local column scaling and gather the results to all ranks
      jmin = minval(jcn(1:nnz))
      jmax = maxval(jcn(1:nnz))
      n_d = jmax - jmin + 1
      allocate(loc_col_scaling(n_d))
      loc_col_scaling(1:n_d) = 1.d-20
      do k=1,nnz
        j = jcn(k) - jmin + 1
        loc_col_scaling(j) = min(max(loc_col_scaling(j),abs(val(k))),1d20)
      enddo          
      do k=1,nnz
        j = jcn(k) - jmin + 1
        val(k) = val(k)/loc_col_scaling(j)
      enddo
      allocate(rcounts(n_cpu_n),displs(n_cpu_n)); rcounts = 0
      rcounts(my_id_n + 1) = n_d
      call MPI_AllReduce(MPI_IN_PLACE, rcounts, n_cpu_n, MPI_INTEGER_ALL, MPI_SUM, comm, ierr)

      displs(1) = 0
      do i = 2, n_cpu_n
        displs(i) = displs(i-1) + rcounts(i-1)
      enddo
      if (allocated(col_scaling)) deallocate(col_scaling)
      allocate(col_scaling(n))
      call MPI_Allgatherv(loc_col_scaling, n_d, MPI_DOUBLE_PRECISION, col_scaling, rcounts, displs, MPI_DOUBLE_PRECISION, comm, ierr)
      deallocate(rcounts,displs,loc_col_scaling)
    
    end subroutine do_column_scaling

#endif
end module mod_pastix

