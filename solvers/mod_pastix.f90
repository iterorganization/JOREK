!> New pastix module to be used with core version
module mod_pastix
#ifdef USE_PASTIX6
  use mod_integer_types
  use iso_c_binding

  type type_PASTIX_SOLVER

    type(C_PTR)                                  :: spm, pastix_data ! sparse solver (distributed)
    type(C_PTR)                                  :: iparm, dparm

    integer                                      :: comm = 0
    real(kind=8), dimension(:), pointer          :: solution_scaling => Null()    !< matrix column scaling to be applied to solution vector
    logical                                      :: initialized = .false.
    logical                                      :: analyzed    = .false.
    logical                                      :: equilibrium = .false.
    logical                                      :: scaled      = .false.
    logical                                      :: refine      = .false.

    integer(kind=int_all), dimension(:), pointer :: loc2glob => null(), glob2loc => null()  ! mapping for column distribution
    real(kind=8), dimension(:), pointer          :: rhs_val => Null()

    integer                                      :: verb = 1                      !< flag for logfile printout (0: no printout)

  end type type_PASTIX_SOLVER

  private
  public :: type_PASTIX_SOLVER, pastix_initialize, pastix_set_mat, pastix_analyze, pastix_factorize, pastix_solve, pastix_finalize

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
      use mod_integer_types
      implicit none

      type(C_PTR), intent(inout) :: spm
      integer :: indx, comm
      integer(kind=int_all) :: n, n_d, nnz, nnz_d, dof
      type(C_PTR), intent(inout) :: rptr, cptr, values, loc2glob, glob2loc
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


  end interface

  contains
!> Initialize PaStiX solver instance
  subroutine pastix_initialize(ptss)
    use, intrinsic :: iso_c_binding
    use mpi_mod
    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    integer ierr

    call ptx_init(ptss%pastix_data,ptss%spm,ptss%iparm,ptss%dparm,ptss%comm)
    call MPI_Barrier(ptss%comm,ierr)

    ptss%initialized = .true.

    return
  end subroutine pastix_initialize

!> Prepare sparse matrix for pastix solver
!! Matrix is distributed column-wise among MPI processes in comm
!! The values are scaled such that the largest value in each column is 1
!! The matrix is converted to CSC format as required by distributed PaStiX
  subroutine pastix_set_mat(ptss, ad_mat, ac_mat, tag)

    use, intrinsic :: iso_c_binding
    use mod_coicsr, only: coicsr, coicsr2
    use sorting_module, only: remove_duplicates, convert2csr
    use mpi_mod
    use data_structure, only: type_SP_MATRIX
    use mod_clock

    implicit none

    type(type_PASTIX_SOLVER)             :: ptss
    type(type_SP_MATRIX)                 :: ad_mat, ac_mat
    integer                              :: tag
    integer                              :: n_cpu, my_id, ierr, comm
    integer(kind=int_all)                :: i, k, nnzg
    integer*8                            :: check_data
    integer                              :: block_size, block_size2, dof
    integer(kind=int_all)                :: nblock, nnz_block
    integer(kind=int_all), allocatable   :: sparskit_work(:)
    integer(kind=int_all)                :: n_d, jmin, jmax
    type(clcktype)                       :: t_itstart, t0, t1, t2, t3
    real*8                               :: tsecond
    type(c_ptr)                          :: irn_c, jcn_c, val_c, loc2glob_c, glob2loc_c
    logical                              :: centralize

    external matrix_split_reduce

    comm = ad_mat%comm

    call MPI_COMM_RANK(comm, my_id, ierr)
    call MPI_COMM_SIZE(comm, n_cpu, ierr)

    if (.not.ptss%equilibrium) then

      call scale_by_cols(ad_mat)
      if (associated(ptss%solution_scaling)) then
        deallocate(ptss%solution_scaling); ptss%solution_scaling => Null()
      endif
      allocate(ptss%solution_scaling(ad_mat%ng))
      do k = 1, ad_mat%ng
        ptss%solution_scaling(k) = ad_mat%column_scaling(k)
      enddo
      ptss%scaled = .true.

    endif

! centralize matrix if it's not distributed by columns; then it will be reduced and distributed
    centralize = (n_cpu.gt.1).and.(.not.ad_mat%col_distributed)
    if (centralize) then
! new ac_mat is allocated, ad_mat is deallocated
      call clck_time(t0)

      call matrix_split_reduce(ad_mat,ac_mat)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (tag .ge. 0)  write(*,FMT_TIMING) tag, '## Elapsed time mpi_gather :', tsecond

    else
! new ac_mat points to ad_mat; no allocation is done
      call ad_mat%move_to(ac_mat, with_data=.true.)
      ac_mat%reduced = .true.
    endif

    call clck_time(t0)

!#ifdef USE_BLOCK
    !block_size = ac_mat%block_size
!#else
    block_size = 1
!#endif
    block_size2 = block_size**2

    nblock   = ac_mat%ng/block_size
    nnz_block = ac_mat%nnz/block_size2

    !if (block_size > 1) then
    !  do i = 1,nnz_block
    !    ac_mat%irn(i) = (ac_mat%irn((i-1)*block_size2+1) - 1)/block_size + 1
    !    ac_mat%jcn(i) = (ac_mat%jcn((i-1)*block_size2+1) - 1)/block_size + 1
    !  enddo
    !endif
    dof = 1

    ! if not already distributed distribute matrix column-wise
    if (ac_mat%col_distributed) then
    ! already column distributed;
      jmin = minval(ac_mat%jcn(1:ac_mat%nnz))
      jmax = maxval(ac_mat%jcn(1:ac_mat%nnz))
      ac_mat%jcn(1:ac_mat%nnz) = ac_mat%jcn(1:ac_mat%nnz) - jmin + ac_mat%indexing
    else
      call distribute_matrix(jmin, jmax, ac_mat)
    endif

    nnzg = 0
    call MPI_Allreduce(ac_mat%nnz,nnzg,1,MPI_INTEGER_ALL,MPI_SUM,comm,ierr)

    n_d = jmax - jmin + 1
    allocate(ptss%loc2glob(n_d))
    do i = 1, n_d
      ptss%loc2glob(i) = i - 1 + jmin;
    enddo

    allocate(ptss%glob2loc(ac_mat%ng)); ptss%glob2loc(1:ac_mat%ng) = 0

    do i = 1, n_d
      ptss%glob2loc(ptss%loc2glob(i) + 1 - ac_mat%indexing)= - my_id - 1;
    enddo

    call MPI_Allreduce(MPI_IN_PLACE, ptss%glob2loc, ac_mat%ng, MPI_INTEGER, MPI_SUM, comm, ierr)

    do i = 1, n_d
      ptss%glob2loc(ptss%loc2glob(i) + 1 - ac_mat%indexing)= ptss%loc2glob(i) - ptss%loc2glob(1)
    enddo

    glob2loc_c = c_loc(ptss%glob2loc); loc2glob_c = c_loc(ptss%loc2glob)

    if (ptss%equilibrium) then
#if (!defined(USEMKL))
      call remove_duplicates(ac_mat%ng, ac_mat%nnz, ac_mat%jcn, ac_mat%irn, ac_mat%val)
#endif
      irn_c = c_loc(ac_mat%irn); jcn_c = c_loc(ac_mat%jcn); val_c = c_loc(ac_mat%val);
      call convert2csr(ac_mat%indexing, n_d, ac_mat%ng, ac_mat%nnz, jcn_c, irn_c, val_c)
    else
      if (allocated(sparskit_work)) deallocate(sparskit_work)
      allocate(sparskit_work(ac_mat%ng + 1))
!#if (defined(USE_BLOCK))
!        call coicsr2(n,nnz_block,val,irn(1:nnz_block),jcn(1:nnz_block),block_size,iwk)
!#else
      call coicsr(ac_mat%ng,ac_mat%nnz,1,ac_mat%val,ac_mat%irn,ac_mat%jcn,sparskit_work)
!#endif
      deallocate(sparskit_work)

      irn_c = c_loc(ac_mat%irn); jcn_c = c_loc(ac_mat%jcn); val_c = c_loc(ac_mat%val);

    endif

    call MPI_Barrier(comm,ierr)

    call ptx_set_mat(ptss%spm, ac_mat%indexing, ac_mat%ng, nnzg, n_d, ac_mat%nnz, dof, irn_c, jcn_c, val_c, loc2glob_c, glob2loc_c, &
                     ptss%comm, ptss%analyzed, ptss%equilibrium)

    return

  end subroutine pastix_set_mat

  subroutine pastix_analyze(ptss)
    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    call ptx_analyze(ptss%pastix_data, ptss%spm)
    return
  end subroutine pastix_analyze

  !> Finalize PaStiX solver instance
  subroutine pastix_finalize(ptss)
    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    call ptx_finalize(ptss%pastix_data, ptss%spm, ptss%iparm, ptss%dparm)

    if (associated(ptss%solution_scaling)) deallocate(ptss%solution_scaling); ptss%solution_scaling => null()
    if (associated(ptss%loc2glob)) deallocate(ptss%loc2glob); ptss%loc2glob => null()
    if (associated(ptss%glob2loc)) deallocate(ptss%glob2loc); ptss%glob2loc => null()

    ptss%comm  = 0

    ptss%scaled      = .false.
    ptss%refine      = .false.
    ptss%equilibrium = .false.

    ptss%initialized = .false.
    ptss%analyzed    = .false.
    ptss%rhs_val     => null() ! not allocated inside pastix module

    return
  end subroutine pastix_finalize

  subroutine pastix_factorize(ptss)
    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    call ptx_factorize(ptss%pastix_data, ptss%spm)

    return
  end subroutine pastix_factorize

  subroutine pastix_solve(ptss,rhs_vec)
    use, intrinsic :: iso_c_binding
    use data_structure, only: type_RHS

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss
    type(type_RHS)                    :: rhs_vec


    type(c_ptr) :: rhsc
    integer(kind=int_all) :: i

    rhsc = c_loc(rhs_vec%val)

    call ptx_solve(ptss%pastix_data, ptss%spm, rhsc, ptss%refine)

    if (ptss%scaled) then
      do i = 1, rhs_vec%n
        rhs_vec%val(i) =  rhs_vec%val(i)/ptss%solution_scaling(i)
      enddo
    endif

    return
  end subroutine pastix_solve

!> Assign local matrix as selected column range from the reduced matrix
  subroutine distribute_matrix(jmin,jmax,ac_mat)
    use data_structure, only: type_SP_MATRIX
    implicit none

    type(type_SP_MATRIX)              :: al_mat, ac_mat

    integer(kind=int_all), intent(out) :: jmin, jmax

    integer :: n_cpu, my_id, comm, ierr
    integer :: i, j, indx

    integer(kind=int_all), allocatable :: dist(:), myelm(:)

    call MPI_Comm_rank(ac_mat%comm, my_id, ierr)
    call MPI_Comm_size(ac_mat%comm, n_cpu, ierr)

    if (allocated(dist)) deallocate(dist)
    allocate(dist(n_cpu + 1))

    ! number of rows/columns per cpu with the last one getting extra
    dist(2:n_cpu + 1) = ac_mat%block_size*((ac_mat%ng/ac_mat%block_size)/n_cpu)
    dist(n_cpu + 1) = dist(n_cpu + 1) + (ac_mat%ng - sum(dist(2:n_cpu + 1)))

    dist(1) = ac_mat%indexing
    do i = 2, n_cpu + 1
      dist(i) = dist(i) + dist(i-1)
    enddo

    jmin = dist(my_id + 1)
    jmax = dist(my_id + 2) - 1

    if (n_cpu.gt.1) then

      allocate(myelm(ac_mat%nnz))
      j = 1
      do i=1, ac_mat%nnz
        if ((ac_mat%jcn(i)>= jmin).and.(ac_mat%jcn(i)<=jmax)) then
          myelm(j) = i
          j = j + 1
        endif
      enddo

      al_mat%nnz = j - 1

      allocate(al_mat%irn(al_mat%nnz),al_mat%jcn(al_mat%nnz),al_mat%val(al_mat%nnz))

      do i = 1, al_mat%nnz
        al_mat%irn(i) = ac_mat%irn(myelm(i))
        al_mat%jcn(i) = ac_mat%jcn(myelm(i)) - dist(my_id + 1) + ac_mat%indexing
        al_mat%val(i) = ac_mat%val(myelm(i))
      enddo

      deallocate(myelm,dist)
      deallocate(ac_mat%irn, ac_mat%jcn, ac_mat%val)

      ac_mat%irn => al_mat%irn
      ac_mat%jcn => al_mat%jcn
      ac_mat%val => al_mat%val
      ac_mat%nnz = al_mat%nnz

    endif

    return

  end subroutine distribute_matrix

#elif USE_PASTIX
#include "pastix_fortran.h"

  use mod_integer_types

  type type_PASTIX_SOLVER
    integer                                      :: comm = 0
    real(kind=8), dimension(:), pointer          :: solution_scaling => Null()    !< matrix column scaling to be applied to solution vector
    logical                                      :: initialized = .false.
    logical                                      :: analyzed    = .false.
    logical                                      :: equilibrium = .false.
    logical                                      :: scaled      = .false.
    logical                                      :: refine      = .false.

    integer(kind=int_all), dimension(:), pointer :: irn => Null()
    integer(kind=int_all), dimension(:), pointer :: jcn => Null()
    real(kind=8), dimension(:), pointer          :: val => Null()
    real(kind=8), dimension(:), pointer          :: rhs_val => Null()

    integer(kind=int_all)                        :: iparm(IPARM_SIZE)
    real*8                                       :: dparm(DPARM_SIZE)

    integer(kind=int_all), dimension(:), pointer :: perm_vars => Null()
    integer(kind=int_all), dimension(:), pointer :: iperm_vars => Null()

    integer(kind=8)                              :: idata    = 0
    integer(kind=int_all)                        :: sym      = API_SYM_NO
    integer(kind=int_all)                        :: iter     = 250
    integer(kind=int_all)                        :: ricar    = 0
    integer(kind=int_all)                        :: iluk     = 3
    integer(kind=int_all)                        :: amalg    = 5
    real*8                                       :: eps      = 1.d-12
    real*8                                       :: pivot    = 1.d-64
    integer(kind=int_all)                        :: maxthrd  = 1024
    integer(kind=int_all)                        :: verb     = API_VERBOSE_NO
    integer(kind=int_all)                        :: facto    = API_FACT_LU
    integer(kind=int_all)                        :: rhs      = 0
    integer(kind=int_all)                        :: nblock   = 0
    integer(kind=8)                              :: block_size = 1
    integer(kind=8)                              :: nthrd = 1

  end type type_PASTIX_SOLVER

  private
  public :: type_PASTIX_SOLVER, pastix_initialize, pastix_set_mat, pastix_analyze, &
            pastix_factorize, pastix_solve, pastix_finalize

  contains

!> Prepares matrix for PaStiX5 solver:
!  rescale columns; convert to csc format
  subroutine pastix_set_mat(ptss, ad_mat, ac_mat, tag)
    use mpi_mod
    use mod_clock
    use mod_coicsr
    use data_structure, only: type_SP_MATRIX
    use tr_module

    implicit none

    type(type_PASTIX_SOLVER)             :: ptss
    type(type_SP_MATRIX)                 :: ad_mat, ac_mat
    integer                              :: tag
    type(clcktype)                       :: t_itstart, t0, t1, t2, t3
    real*8                               :: tsecond
    integer                              :: n_cpu, my_id, ierr, comm
    integer                              :: i, j
    integer(kind=int_all)                :: k, nnz
    integer*8                            :: check_data
    integer(kind=int_all)                :: nblock, nnz_block, block_size, block_size2
    integer(kind=int_all), allocatable   :: sparskit_work(:)
    logical                              :: distributed

    external matrix_split_reduce, matrix_split_bcast

    comm = ad_mat%comm

    call MPI_COMM_RANK(comm, my_id, ierr)
    call MPI_COMM_SIZE(comm, n_cpu, ierr)

    if (.not.ptss%equilibrium) then

      call scale_by_cols(ad_mat)
      if (associated(ptss%solution_scaling)) then
        deallocate(ptss%solution_scaling); ptss%solution_scaling => Null()
      endif
      allocate(ptss%solution_scaling(ad_mat%ng))
      do k = 1, ad_mat%ng
        ptss%solution_scaling(k) = ad_mat%column_scaling(k)
      enddo
      ptss%scaled = .true.

    endif

! centralize matrix if it's not distributed by columns; then it will be reduced and distributed
! use reduce if not row_distributed (meaning matrix is on my_id_n==0)
    distributed = ad_mat%row_distributed.or.ad_mat%col_distributed
    if ((n_cpu.gt.1).and.(.not.distributed)) then
    ! matrix comes from communication, thus must be broadcasted
      call clck_time(t0)

      call matrix_split_bcast(ad_mat,ac_mat)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (tag .ge. 0)  write(*,FMT_TIMING) tag, '## Elapsed time mpi_bcast :', tsecond


    elseif ((n_cpu.gt.1).and.(distributed)) then
    ! matrix comes from direct construction, thus must be all-gathered
      call clck_time(t0)

      call matrix_split_reduce(ad_mat,ac_mat)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (tag .ge. 0)  write(*,FMT_TIMING) tag, '## Elapsed time mpi_gather :', tsecond

    else
    ! matrix is already reduced (n_cpu = 1)
      call ad_mat%move_to(ac_mat, with_data=.true.)
      ac_mat%reduced = .true.
    endif

    call clck_time(t0)

#ifdef USE_BLOCK
    block_size = ac_mat%block_size
#else
    block_size = 1
#endif
    block_size2 = block_size**2
    nblock   = ac_mat%ng/block_size
    nnz_block = ac_mat%nnz/block_size2

    if (block_size > 1) then
      do i = 1,nnz_block
        ac_mat%irn(i) = (ac_mat%irn((i-1)*block_size2+1) - 1)/ac_mat%block_size + 1
        ac_mat%jcn(i) = (ac_mat%jcn((i-1)*block_size2+1) - 1)/ac_mat%block_size + 1
      enddo
    endif

    allocate(sparskit_work(nblock+1))

    call coicsr2(nblock,nnz_block,ac_mat%val,ac_mat%irn(1:nnz_block),ac_mat%jcn(1:nnz_block),block_size,sparskit_work)

    deallocate(sparskit_work)

    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (tag .ge. 0)  write(*,FMT_TIMING) tag, '## Elapsed time coicsr :', tsecond

    ! End of matrix preparation
    ptss%irn => ac_mat%irn
    ptss%jcn => ac_mat%jcn
    ptss%val => ac_mat%val
    ptss%nblock = nblock
    ptss%block_size = block_size

    if (ptss%equilibrium) then
    ! combine duplicated values
      nnz = ac_mat%jcn(ac_mat%ng + 1) - 1
      call pastix_fortran_checkmatrix(check_data, ac_mat%comm, &
       Int1, ptss%sym, Int1, ac_mat%ng, ac_mat%jcn, ac_mat%irn, ac_mat%val, -Int1, Int1)

      ac_mat%nnz = ac_mat%jcn(ac_mat%ng+1) - 1
      if (ac_mat%nnz < nnz) then
         call pastix_fortran_checkmatrix_end(check_data, Int1, ac_mat%irn, ac_mat%val, Int1)
      endif
    endif

    return

  end subroutine pastix_set_mat

!> Initializes PaStiX5 solver
  subroutine pastix_initialize(ptss)
    use mpi_mod
    use mod_clock
    use phys_module, only: pastix_maxthrd

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss
    integer                           :: my_id, ierr
    type(clcktype)                    :: t_itstart, t0, t1, t2, t3
    real*8                            :: tsecond

    ptss%maxthrd = pastix_maxthrd
    call MPI_COMM_RANK(ptss%comm, my_id, ierr)

    call pastix_init_nthreads(ptss)

    if (my_id .eq. 0) write(*,'(A,i5)') ' PastiX n_threads : ', ptss%nthrd

    if (associated(ptss%perm_vars)) deallocate(ptss%perm_vars); ptss%perm_vars  => Null()
    if (associated(ptss%iperm_vars)) deallocate(ptss%iperm_vars); ptss%iperm_vars => Null()

    allocate(ptss%perm_vars(ptss%nblock))
    allocate(ptss%iperm_vars(ptss%nblock))
    ptss%perm_vars(1:ptss%nblock) = 0
    ptss%iperm_vars(1:ptss%nblock) = 0

    ptss%iparm(IPARM_MODIFY_PARAMETER)  = API_NO          ! insert default values
    ptss%iparm(IPARM_START_TASK)        = API_TASK_INIT   ! initializse
    ptss%iparm(IPARM_END_TASK)          = API_TASK_INIT

    if (my_id .eq. 0) write(*,*) "Initializing PaStiX solver"

    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)

    ptss%iparm(IPARM_VERBOSE)               = ptss%verb
    ptss%iparm(IPARM_ITERMAX)               = ptss%iter

    ptss%iparm(IPARM_FACTORIZATION)         = ptss%facto
    ptss%iparm(IPARM_INCOMPLETE)            = ptss%ricar
    ptss%iparm(IPARM_LEVEL_OF_FILL)         = ptss%iluk
    ptss%dparm(DPARM_EPSILON_REFINEMENT)    = ptss%eps
    ptss%dparm(DPARM_EPSILON_MAGN_CTRL)     = ptss%pivot
    ptss%iparm(IPARM_RHS_MAKING)            = ptss%rhs
    ptss%iparm(IPARM_SYM)                   = ptss%sym
    ptss%iparm(IPARM_AMALGAMATION_LEVEL)    = ptss%amalg

#ifdef FUNNELED
    ptss%iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_FUNNELED
#else
    ptss%iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_MULTIPLE
#endif

    ptss%initialized = .true.

    return

  end subroutine pastix_initialize

  subroutine pastix_finalize(ptss)
    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    ptss%iparm(IPARM_START_TASK) = API_TASK_CLEAN
    ptss%iparm(IPARM_END_TASK)   = API_TASK_CLEAN

    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)

    if (associated(ptss%perm_vars)) deallocate(ptss%perm_vars); ptss%perm_vars  => Null()
    if (associated(ptss%iperm_vars)) deallocate(ptss%iperm_vars); ptss%iperm_vars => Null()

    ptss%irn => Null()
    ptss%jcn => Null()
    ptss%val => Null()
    ptss%rhs_val => Null()

    ptss%idata = 0
    ptss%comm  = 0

    ptss%scaled      = .false.
    ptss%refine      = .false.
    ptss%equilibrium = .false.

    ptss%initialized = .false.
    ptss%analyzed    = .false.


    return

  end subroutine pastix_finalize

!> Performs matrix analyzis/reordering with PaStiX5 solver
  subroutine pastix_analyze(ptss)
    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    ptss%iparm(IPARM_DOF_NBR) = ptss%block_size
    ptss%iparm(IPARM_THREAD_NBR) = ptss%nthrd

    ptss%iparm(IPARM_START_TASK) = API_TASK_ORDERING
    ptss%iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)

    ptss%analyzed = .true.

    return

  end subroutine pastix_analyze

!> Performs matrix LU factorization with PaStiX5 solver
  subroutine pastix_factorize(ptss)

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss

    ptss%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
    ptss%iparm(IPARM_END_TASK)   = API_TASK_NUMFACT

    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, ptss%rhs_val, int1, ptss%iparm, ptss%dparm)

    return

  end subroutine pastix_factorize

!> Finds the solution with PaStiX5 solver for given RHS
  subroutine pastix_solve(ptss,rhs_vec)
    use data_structure, only: type_RHS

    implicit none

    type(type_PASTIX_SOLVER)          :: ptss
    type(type_RHS)                    :: rhs_vec

    integer(kind=int_all)             :: i
    ! dummy arguments to be used with pastix interface
    integer(kind=int_all), pointer    :: irn(:), jcn(:)
    real(kind=8), pointer             :: val(:)

    ptss%iparm(IPARM_START_TASK) = API_TASK_SOLVE
    ptss%iparm(IPARM_END_TASK)   = API_TASK_SOLVE
    if (ptss%refine) ptss%iparm(IPARM_END_TASK) = API_TASK_REFINE

    call pastix_fortran(ptss%idata, ptss%comm, ptss%nblock, ptss%jcn, ptss%irn, ptss%val, &
                        ptss%perm_vars, ptss%iperm_vars, rhs_vec%val, int1, ptss%iparm, ptss%dparm)

    if (ptss%scaled) then
      do i=1,rhs_vec%n
        rhs_vec%val(i) =  rhs_vec%val(i)/ptss%solution_scaling(i)
      enddo
    endif

    return

  end subroutine pastix_solve


  subroutine pastix_init_nthreads(ptss)
    use mpi_mod
#ifdef _OPENMP
    use omp_lib
#endif

    implicit none
    type(type_PASTIX_SOLVER) :: ptss
    integer :: nbthreads

#ifdef _OPENMP
    !$OMP PARALLEL shared(nbthreads)
    !$OMP master
    call omp_set_dynamic(.false.)
    nbthreads = omp_get_num_threads()
    !$OMP end master
    !$OMP barrier
    !$OMP end PARALLEL
#else
    nbthreads = 1
#endif

    ptss%nthrd = max(min(ptss%maxthrd,nbthreads),1)

    return
  end subroutine pastix_init_nthreads

#endif
end module mod_pastix

