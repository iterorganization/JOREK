!> New strumpack module to be used with core version
module mod_strumpack
#ifdef USE_STRUMPACK
  use iso_c_binding
  use mod_integer_types
  use mpi_mod

  implicit none

  type type_STRUMPACK_SOLVER
    type(c_ptr) :: sscp ! STRUMPACK sparse solver c pointer
    !type(c_ptr) :: distr
    integer(kind=C_INT_ALL), pointer :: distr(:)
    logical                          :: initialized = .false.
    logical                          :: analyzed = .false.
    logical                          :: equilibrium = .false.
    integer                          :: comm = 0
  end type type_STRUMPACK_SOLVER

  private
  public :: type_STRUMPACK_SOLVER, &
            strumpack_init, strumpack_set_mat, strumpack_analyze, &
            strumpack_factorize, strumpack_solve, strumpack_finalize

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
      type(c_ptr) :: irn, jcn, val
      type(c_ptr) :: dist
      !integer(kind=C_INT_ALL), dimension(:), pointer, intent(in) :: dist
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
      type(c_ptr), intent(inout) :: spss, rhs
      type(c_ptr) :: dist
      !integer(kind=C_INT_ALL), dimension(:), pointer, intent(in) :: dist
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
    
    subroutine strumpack_init(spss,comm)
      implicit none

      integer comm,ierr, n_cpu
      type(type_STRUMPACK_SOLVER) spss
      
      spss%comm = comm

      call spk_init(spss%sscp,spss%comm)
      
      call MPI_COMM_SIZE(spss%comm, n_cpu, ierr)
      allocate(spss%distr(n_cpu+1))
      
      call MPI_Barrier(spss%comm,ierr)

      return
    end subroutine strumpack_init  
    
    
    subroutine strumpack_set_mat(spss,a_mat)
      use, intrinsic :: iso_c_binding
      use sorting_module, only : remove_duplicates, convert2csr, convert_sorting
      use mod_integer_types
      use data_structure, only: type_SP_MATRIX

      implicit none
      
      type(type_STRUMPACK_SOLVER)       :: spss
      type(type_SP_MATRIX)              :: a_mat

      integer ierr
      integer(kind=C_INT_ALL), dimension(:), pointer :: irn_d, jcn_d
      real(kind=C_DOUBLE),  dimension(:), pointer :: val_d

      integer(kind=C_INT_ALL), allocatable, target :: distr(:)
      
      integer :: rank, n_cpu
      integer(kind=int_all) :: nnz_d, n_d, i, j, imin, imax

      integer(kind=int_all), dimension(:), pointer :: myelm
      logical :: upd, dflag, eql
      type(c_ptr) :: irn_c, jcn_c, val_c, dist_c
      
      upd = spss%analyzed
      dflag = a_mat%row_distributed
      eql = spss%equilibrium

      call MPI_COMM_RANK(spss%comm, rank, ierr)
      call MPI_COMM_SIZE(spss%comm, n_cpu, ierr)
      
      !write(*,*) "n_cpu, a_mat%ng, a_mat%nnz", n_cpu, a_mat%ng, a_mat%nnz
      
      allocate(distr(n_cpu+1))

      if ((.not. dflag).and.(n_cpu.gt.1)) then
        ! distribute rows between n_cpu
        call distribute_rows(a_mat,n_cpu,distr)
        if (rank.eq.0) write(*,*) "Matrix is not row-distributed. Distributing now."

        allocate(myelm(a_mat%nnz))
        j = 1
        do i=1, a_mat%nnz
          if ((a_mat%irn(i) >= distr(rank+1)).and.(a_mat%irn(i) <= (distr(rank+2)-1))) then
            myelm(j) = i
            j = j + 1
          endif
        enddo

        nnz_d = j - 1
        n_d = distr(rank+2) - distr(rank+1) ! number of local rows

        allocate(irn_d(nnz_d), jcn_d(nnz_d), val_d(nnz_d))

        do i = 1, nnz_d
          irn_d(i) = a_mat%irn(myelm(i)) - distr(rank+1) + a_mat%indexing    ! irn starts from index
          jcn_d(i) = a_mat%jcn(myelm(i))                          ! jcn remains the same
          val_d(i) = a_mat%val(myelm(i))
        enddo

        deallocate(a_mat%irn,a_mat%jcn,a_mat%val)
        
        a_mat%irn => irn_d
        a_mat%jcn => jcn_d
        a_mat%val => val_d
        a_mat%row_distributed = .true.

        deallocate(myelm)

      elseif (dflag.and.(n_cpu.gt.1)) then
        ! get row distribution from irn in case of pre-distributed matrix
        if (allocated(distr)) deallocate(distr)
        allocate(distr(n_cpu+1))
        
        distr(1:n_cpu+1) = 0
        imin = minval(a_mat%irn(1:a_mat%nnz))
        imax = maxval(a_mat%irn(1:a_mat%nnz))
                
        distr(rank+1) = imin

        if (rank.eq.(n_cpu-1)) distr(rank+2) = imax + 1
        call MPI_Allreduce(MPI_IN_PLACE,distr,n_cpu+1,MPI_INTEGER,MPI_SUM,spss%comm,ierr)

        ! check for consistency
        ierr = 0
        if ((distr(1).ne.1)) ierr = 1
        do i = 2, n_cpu+1
          if (.not.(distr(i)>distr(i-1))) ierr = 1
        enddo

        if (ierr.ne.0) then
          write(*,*) "Error in harmonic matrix distribution"
          call exit(2)
        endif

        n_d = distr(rank+2) - distr(rank+1)
        nnz_d = a_mat%nnz
        a_mat%irn(1:nnz_d) = a_mat%irn(1:nnz_d) - imin + a_mat%indexing ! irn starts with indx
        
      elseif (n_cpu.eq.1) then
        call distribute_rows(a_mat,1,distr)
        n_d = a_mat%ng
        nnz_d = a_mat%nnz
        write(*,*) "n_d, nnz_d", n_d, nnz_d

      endif
      
      if (eql) then
        call remove_duplicates(a_mat%ng,a_mat%nnz,a_mat%irn,a_mat%jcn,a_mat%val)
        n_d = a_mat%ng
        nnz_d = a_mat%nnz
      endif
      
      if (a_mat%indexing.eq.1) then
        a_mat%irn(1:nnz_d) = a_mat%irn(1:nnz_d) - a_mat%indexing;
        a_mat%jcn(1:nnz_d) = a_mat%jcn(1:nnz_d) - a_mat%indexing;
        distr(1:n_cpu+1) = distr(1:n_cpu+1) - a_mat%indexing
        a_mat%indexing = 0
      endif
      
#if (!defined(USEMKL))
      call convert_sorting(nnz_d, a_mat%irn, a_mat%jcn, a_mat%val, a_mat%block_size, a_mat%indexing)
#endif

      irn_c = c_loc(a_mat%irn); jcn_c = c_loc(a_mat%jcn); val_c = c_loc(a_mat%val); dist_c = c_loc(distr)
      
#if (defined(USEMKL))
      call convert2csr(a_mat%indexing, n_d, a_mat%ng, nnz_d, irn_c, jcn_c, val_c)
#endif

      !spss%distr = c_loc(distr)
      
      spss%distr(1:n_cpu+1) = distr(1:n_cpu+1)
      deallocate(distr)
      dist_c = c_loc(spss%distr)
      
      call spk_set_mat(n_d,dist_c,irn_c,jcn_c,val_c,spss%sscp,spss%comm,upd)

      call MPI_Barrier(spss%comm,ierr)


      return
    end subroutine strumpack_set_mat
    
    subroutine strumpack_analyze(spss)
      use data_structure, only: type_SP_MATRIX

      implicit none
      
      type(type_STRUMPACK_SOLVER)       :: spss
      integer ierr

      call spk_reord(spss%sscp,spss%comm)
      call MPI_Barrier(spss%comm,ierr)
      spss%analyzed = .true.

      return
    end subroutine strumpack_analyze 
   
    subroutine strumpack_factorize(spss)
      use data_structure, only: type_SP_MATRIX

      implicit none
      
      type(type_STRUMPACK_SOLVER)       :: spss
      integer ierr

      call spk_fact(spss%sscp,spss%comm)
      call MPI_Barrier(spss%comm,ierr)

      return
    end subroutine strumpack_factorize
   
    subroutine strumpack_solve(spss, rhs_vec)
      use data_structure, only: type_SP_MATRIX, type_RHS
      use, intrinsic :: iso_c_binding

      implicit none
      
      type(type_STRUMPACK_SOLVER)   :: spss
      type(type_RHS)                :: rhs_vec
      type(c_ptr)                   :: rhs_c
      type(c_ptr)                   :: dist_c
      integer :: ierr
      integer(kind=C_INT_ALL), allocatable, target :: dist(:)
      integer :: i, n_cpu
      
      call MPI_COMM_SIZE(spss%comm, n_cpu, ierr) 

      rhs_c = c_loc(rhs_vec%val);
      dist_c = c_loc(spss%distr)

      call spk_solve(rhs_vec%n, dist_c, rhs_c, spss%sscp, spss%comm)
      call MPI_Barrier(spss%comm,ierr)

      return
    end subroutine strumpack_solve

!> Finalize strumpack solver instance
    subroutine strumpack_finalize(spss)
      implicit none

      type(type_STRUMPACK_SOLVER)   :: spss

      if (associated(spss%distr)) then
        deallocate(spss%distr); spss%distr => null()
      endif

      call spk_finalize(spss%sscp, spss%comm)

      spss%initialized = .false.
      spss%analyzed    = .false.
      spss%equilibrium = .false.
      spss%comm        = 0

      return
    end subroutine strumpack_finalize
   
!> Distribute rows between members of MPI group   
    subroutine distribute_rows(a_mat,n_cpu,distr)
      
      use, intrinsic :: iso_c_binding
      use data_structure, only: type_SP_MATRIX
      use mod_integer_types
      implicit none

      type(type_SP_MATRIX)                 :: a_mat
      integer, intent(in)                  :: n_cpu
      
      integer(kind=C_INT_ALL), allocatable :: nr(:)
      integer(kind=C_INT_ALL), allocatable :: distr(:)
      integer :: ierr, i
    
      
      allocate(nr(n_cpu))

      nr = a_mat%block_size*((a_mat%ng/a_mat%block_size)/n_cpu)
      nr(n_cpu) = nr(n_cpu) + (a_mat%ng - sum(nr))

      distr(1) = 1
      do i=1, n_cpu
        distr(i+1)= distr(i) + nr(i)
      enddo

      deallocate(nr)

      return
    end subroutine distribute_rows

#endif
end module mod_strumpack

