! module dedicated to memory tracing and
! debug tracing
module tr_module

  use mod_integer_types

  implicit none
  interface tr_register_mem
     module procedure tr_register_mem_int4, tr_register_mem_int8
  end interface

  interface tr_unregister_mem
     module procedure tr_unregister_mem_int4, tr_unregister_mem_int8
  end interface

  !*** surdefinition for allocation ***
  interface tr_allocate
     module procedure &
          tr_allocate1d_ch,                 &
          tr_allocate1d_i, tr_allocate1d_d, &
          tr_allocate1d_c, tr_allocate1d_ll,&
          tr_allocate1d_ld,tr_allocate1d_il,&
          tr_allocate1d_li,                 &
          tr_allocate2d_i, tr_allocate2d_d, &
          tr_allocate2d_c, tr_allocate2d_ll,&
          tr_allocate3d_i, tr_allocate3d_d, &
          tr_allocate3d_c,                  &
          tr_allocate4d_d
  end interface

  !*** surdefinition for deallocation ***
  interface tr_deallocate
     module procedure &
          tr_deallocate1d_ch,                   &
          tr_deallocate1d_i, tr_deallocate1d_d, &
          tr_deallocate1d_c, tr_deallocate1d_l, &
          tr_deallocate2d_i, tr_deallocate2d_d, &
          tr_deallocate2d_c, tr_deallocate2d_l, &
          tr_deallocate3d_i, tr_deallocate3d_d, &
          tr_deallocate3d_c,                    &
          tr_deallocate4d_d
  end interface

  !*** surdefinition for allocation ***
  interface tr_allocatep
     module procedure tr_allocatep1d_i,       &
          tr_allocatep1d_d, tr_allocatep1d_ll,&
          tr_allocatep1d_ld,tr_allocatep1d_il,&
          tr_allocatep1d_li,                  &
          tr_allocatep2d_d, tr_allocatep2d_i, &
          tr_allocatep2d_ll,                  &
          tr_allocatep3d_d, tr_allocatep3d_i, &
          tr_allocatep4d_d, tr_allocatep4d_i
  end interface

  !*** surdefinition for deallocation ***
  interface tr_deallocatep
     module procedure tr_deallocatep1d_i, tr_deallocatep1d_l, &
          tr_deallocatep1d_d, tr_deallocatep2d_i, &
          tr_deallocatep2d_l,                     &
          tr_deallocatep3d_d, tr_deallocatep2d_d, &
          tr_deallocatep3d_i, tr_deallocatep4d_d, &
          tr_deallocatep4d_i
  end interface

  !*** surdefinition for deallocation ***
  interface tr_debug_write
     module procedure tr_debug_writes, tr_debug_writei, tr_debug_writel, tr_debug_writef
  end interface

  real*8 :: precond_mem = 0
  ! size of real numbers used in jorek
  integer, public, parameter :: RKIND = 8
  ! processor identity
  integer, public :: gmy_id
  ! nb of processors
  integer, private :: nbprocs
  character(LEN=13), private :: &
       trace_file = "trace    .out"
  integer, private :: uout_mem = 30

  real(RKIND) :: myreal
  integer     :: myint
  complex     :: mycomp
  public 
  integer, PARAMETER :: CAT_UNKNOWN = 0
  integer, PARAMETER :: CAT_MATELEM = 1
  integer, PARAMETER :: CAT_GRID    = 2
  integer, PARAMETER :: CAT_FEM     = 3
  integer, PARAMETER :: CAT_DMATRIX = 4
  integer, PARAMETER :: CAT_GLOBMAT = 5
  integer, PARAMETER :: CAT_PRECOND = 6
  integer, PARAMETER :: CAT_GMRES   = 7
  integer, PARAMETER :: MIN_CAT     = CAT_UNKNOWN
  integer, PARAMETER :: MAX_CAT     = CAT_GMRES
  CHARACTER(LEN=10) :: &
    cat_name(MIN_CAT:MAX_CAT) = (/ "UNKNOWN   ", "MATELEM   ", "GRID      ", &
    "FEM       ",  "DMATRIX   ", "GLOBMAT   ", "PRECOND   ", "GMRES     " /)
  ! used for memory size calculation
  integer*8, private :: max_allocate
  integer*8, private :: nb_allocate(MIN_CAT:MAX_CAT)

  !******************************
contains
  !******************************

  subroutine tr_vdump(filename,mat,nnz)
    REAL*8, dimension(:) :: mat
    character(len=*) :: filename
    character(len=len(filename)+10) :: filename2
    INTEGER(kind=int_all) :: nnz
    
#ifdef NORMTRACE
    INTEGER i, j
    LOGICAL file_exists
    INQUIRE(FILE=filename, EXIST=file_exists)
    if (file_exists) then
       j = 0
       do while (file_exists)
          write(filename2, '(A,A1,I6.6)') filename, "_", j
          INQUIRE(FILE=filename2, EXIST=file_exists)
          j = j + 1
       end do
       CALL RENAME(filename, filename2)
    end if
    open(11, file = filename, status = 'REPLACE', form = 'FORMATTED')
    print *, "mat(1)", mat(1)
    write (*,"(A6,E20.12)") "mat(1)", mat(1)
    do i = 1, nnz
       write (11,"(I10,E20.12)") i, mat(i)
    end do
    print *, "mat(1)", mat(1)
    close(11)
#endif
  end subroutine tr_vdump


  subroutine tr_vnorms(prefix,mat,nnz)
    use mpi_mod
    REAL*8, dimension(:) :: mat
    character(len=*) :: prefix
    INTEGER(kind=int_all) :: nnz

#ifdef NORMTRACE
    REAL*8  :: lnorms(1:3)
    REAL*8  :: l1, l2, linf, absv
    character(len=92) :: bufstring
    INTEGER :: ierr
    INTEGER :: i
    l1 = 0._8
    l2 = 0._8
    linf = 0._8
    do i = 1, nnz
       absv = abs(mat(i))
       l2   = l2   + mat(i)*mat(i) 
       l1   = l1   + absv
       if (absv .gt. linf) linf = absv
    end do
    call MPI_AllReduce(l1,lnorms(1),1,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(l2,lnorms(2),1,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(linf,lnorms(3),1,MPI_DOUBLE_PRECISION,&
         MPI_MAX,MPI_COMM_WORLD,ierr)
    write(bufstring,'(A20,A4,1X,3E22.13)') trim(prefix), ":gn:", lnorms(1:3)
    call tr_write(bufstring)
#endif
  end subroutine tr_vnorms

  subroutine tr_locvnorms(prefix,mat,nnz)
    use mpi_mod
    REAL*8, dimension(:) :: mat
    character(len=*) :: prefix
    INTEGER(kind=int_all) :: nnz
#ifdef NORMTRACE
    REAL*8  :: lnorms(1:3)
    REAL*8  :: l1, l2, linf, absv
    character(len=92) :: bufstring
    INTEGER :: ierr
    INTEGER :: i
    l1 = 0._8
    l2 = 0._8
    linf = 0._8
    do i = 1, nnz
       absv = abs(mat(i))
       l2   = l2   + mat(i)*mat(i) 
       l1   = l1   + absv
       if (absv .gt. linf) linf = absv
    end do
    lnorms(1:3) = (/ l1, l2, linf /)
    write(bufstring,'(A20,A4,1X,3E22.13)') trim(prefix), ":ln:", lnorms(1:3)
    call tr_write(bufstring)
#endif
  end subroutine tr_locvnorms

  subroutine tr_locvnorms_cmplx(prefix,mat,nnz)
    use mpi_mod
    double complex, dimension(:) :: mat
    character(len=*) :: prefix
    INTEGER :: nnz
#ifdef NORMTRACE
    double complex  :: lnorms(1:3)
    double complex  :: l1, l2, linf, absv
    character(len=92) :: bufstring
    INTEGER :: ierr
    INTEGER :: i

    l1 = COMPLEX(0._8,0._8)
    l2 = COMPLEX(0._8,0._8)
    linf = COMPLEX(0._8,0._8)
    do i = 1, nnz
       absv = abs(mat(i))
       l2   = l2   + mat(i)*mat(i) 
       l1   = l1   + absv
       if (absv .gt. linf) linf = absv
    end do
    lnorms(1:3) = (/ l1, l2, linf /)
    write(bufstring,'(A20,A4,1X,3E22.13)') trim(prefix), ":ln:", lnorms(1:3)
    call tr_write(bufstring)
#endif
  end subroutine tr_locvnorms_cmplx



  !---------------------------------------- 
  ! Init target file in each processor
  !----------------------------------------
  subroutine tr_meminit(pmy_id, pnbprocs)
    integer, intent(in) :: pmy_id, pnbprocs
    logical fexist
    gmy_id = pmy_id
    nbprocs = pnbprocs
    write(trace_file(6:9),'(I4.4)') gmy_id
    inquire(file=trace_file,exist=fexist)
    if (fexist) then
       open(uout_mem, file = trace_file, status = 'OLD', position = 'APPEND', &
            form = 'FORMATTED')
    else
       open(uout_mem, file = trace_file, status = 'NEW', &
            form = 'FORMATTED')
    end if
    write(uout_mem,*) '### meminit ### '
    call flush_it(uout_mem)
    close(uout_mem)
    nb_allocate(:) = 0
  end subroutine tr_meminit

  !---------------------------------------- 
  ! Reset file
  !----------------------------------------
  subroutine tr_resetfile()
    open(uout_mem, file = trace_file, status = 'REPLACE', &
         form = 'FORMATTED')
    write(uout_mem,*) '### meminit ### '
    call flush_it(uout_mem)
    close(uout_mem)
  end subroutine tr_resetfile


  !---------------------------------------- 
  ! Write special string in file trace_file
  !----------------------------------------
  subroutine tr_write(string)
    character*(*) string
    open(uout_mem, file = trace_file, status = 'OLD', &
         position = 'APPEND', form = 'FORMATTED')
    write(uout_mem,'(A)') string
    call flush_it(uout_mem)
    close(uout_mem)
  end subroutine tr_write

  !---------------------------------------- 
  ! Write debug remark in file trace_file
  !----------------------------------------
  subroutine tr_debug_writes(string)
    character*(*)           :: string
    character(len=1024)     :: bufstring
    write(bufstring,'(A)')string
    call tr_write("### "//trim(adjustl(bufstring))//" ###")
  end subroutine tr_debug_writes


  !---------------------------------------- 
  ! Write debug remark in file trace_file
  !----------------------------------------
  subroutine tr_debug_writei(string, int_var)
    character*(*)           :: string
    integer                 :: int_var
    character(len=1024)     :: bufstring
    write(bufstring,'(A,I20)')string,int_var
    call tr_write("### "//trim(adjustl(bufstring))//" ###")
  end subroutine tr_debug_writei

  !---------------------------------------- 
  ! Write debug remark in file trace_file
  !----------------------------------------
  subroutine tr_debug_writel(string, int_var)
    character*(*)           :: string
    integer*8               :: int_var
    character(len=1024)     :: bufstring
    write(bufstring,'(A,I20)')string,int_var
    call tr_write("### "//trim(adjustl(bufstring))//" ###")
  end subroutine tr_debug_writel

  !---------------------------------------- 
  ! Write debug remark in file trace_file
  !----------------------------------------
  subroutine tr_debug_writef(string, float_var)
    character*(*)           :: string
    real*8                  :: float_var
    character(len=1024)     :: bufstring
    write(bufstring,'(A,E20.11)')string,float_var
    call tr_write("### "//trim(adjustl(bufstring))//" ###")
  end subroutine tr_debug_writef


  !-------------------------------------------
  ! Write memory in the file trace_file (allocate)
  !-------------------------------------------
  subroutine tr_memwriteadd(size_array,type_name,var_name,category)
    integer*8    , intent(in)  :: size_array
    character*(*), intent(in)  :: type_name
    character*(*), intent(in)  :: var_name
    integer      , intent(in)  :: category
#ifdef MEMTRACE
    open(uout_mem, file = trace_file, status = 'OLD', &
         position = 'APPEND', form = 'FORMATTED')
    write(uout_mem,'(A10,I15,A3,A15,1X,A15,A20,5X,I20)') &
         'Add', &
         size_array,' B ',type_name,cat_name(category),&
         var_name,SUM(nb_allocate(MIN_CAT:MAX_CAT))
    close(uout_mem)    
#endif
  end subroutine tr_memwriteadd

  !-------------------------------------------
  ! Write memory in the file trace_file (deallocate)
  !-------------------------------------------
  subroutine tr_memwritedel(var_name,category)
    character*(*), intent(in) :: var_name
    integer      , intent(in) :: category

#ifdef MEMTRACE
    open(uout_mem, file = trace_file, status = 'OLD', &
         position = 'APPEND', form = 'FORMATTED')
    write(uout_mem,'(30X,A20,1X,A10,A5,I20)') var_name,cat_name(category),' Supp', &
      SUM(nb_allocate(MIN_CAT:MAX_CAT))
    close(uout_mem)    
#endif
  end subroutine tr_memwritedel


  !-------------------------------------------
  ! Get memory usage in looking at /proc/self/status
  ! Return an integer representing the consumption in KB
  !-------------------------------------------
  integer*8 function get_memory_inkb(string_pattern)
    character(len=*), intent(in) :: string_pattern
    character(len=1024) :: buffer, string, adj
    integer :: pos, ierr
    integer, parameter :: fh = 15
    integer :: ios 
    integer :: line 
    ios = 0
    line = 0
    get_memory_inkb = 0
    open(fh, file='/proc/self/status', iostat=ierr)
    if (ierr .ne. 0) get_memory_inkb = -1
    ! ios is negative if an end of record condition is encountered or if                                                                                                      
    ! an endfile condition was detected.  It is positive if an error was                                                                                                      
    ! detected.  ios is zero otherwise.                                                                                                                                       
    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1
          ! Find the first instance of whitespace.                                                                                                                            
          ! Split first field and end of line                                                                                                                                 
          pos = scan(buffer,':')
          string = buffer(1:pos-1)
          if (trim(adjustl(string)) .eq. trim(adjustl(string_pattern))) then
             buffer = buffer(pos+1:)
             pos = scan(buffer,'0123456789')
             adj = buffer(pos:)
             pos = scan(adj,'k')
             string = adj(1:pos-2)
             read(string,'(I20)')get_memory_inkb
          end if
       end if
    end do
    close(fh)
    return 
  end function get_memory_inkb

  subroutine tr_set_precondmem(memused)
    real*8 memused
    precond_mem = memused
  end subroutine tr_set_precondmem


  subroutine tr_register_mem_int4(mem_in_bytes,var_name,ocategory)
    integer*4          , intent(in) :: mem_in_bytes
    character*(*)      , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    call tr_memwriteadd(int(mem_in_bytes,8),'unknown type',var_name,category)
    nb_allocate(category) = nb_allocate(category) + mem_in_bytes
  end subroutine tr_register_mem_int4

  subroutine tr_register_mem_int8(mem_in_bytes,var_name,ocategory)
    integer*8          , intent(in) :: mem_in_bytes
    character*(*)      , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    call tr_memwriteadd(mem_in_bytes,'unknown type',var_name,category)
    nb_allocate(category) = nb_allocate(category) + mem_in_bytes
  end subroutine tr_register_mem_int8

  subroutine tr_unregister_mem_int4(mem_in_bytes,var_name,ocategory)
    integer*4          , intent(in) :: mem_in_bytes
    character*(*)      , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    nb_allocate(category) = nb_allocate(category) - mem_in_bytes
    call tr_memwritedel(var_name,category)
  end subroutine tr_unregister_mem_int4

  subroutine tr_unregister_mem_int8(mem_in_bytes,var_name,ocategory)
    integer*8          , intent(in) :: mem_in_bytes
    character*(*)      , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    nb_allocate(category) = nb_allocate(category) - mem_in_bytes
    call tr_memwritedel(var_name,category)
  end subroutine tr_unregister_mem_int8

  !---------------------------------------- 
  ! memory allocation for a 1D array
  !----------------------------------------
  subroutine tr_allocatep1d_i(array1d,begin_dim1,end_dim1,var_name,ocategory)
    integer, dimension(:)    , pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep1d_i

  subroutine tr_allocatep1d_d(array1d,begin_dim1,end_dim1,var_name,ocategory,pzeroing)
    real(RKIND), dimension(:), pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    logical        , optional, intent(in) :: pzeroing
    integer        , optional, intent(in) :: ocategory
    logical :: zeroing
    integer :: category

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (.not. present(pzeroing)) then
      zeroing = .true.
    else
      zeroing = pzeroing
    end if
    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array1D',var_name,category)
    if (err.eq.0) then
       if (zeroing) then
          do i1 = begin_dim1,end_dim1
             array1d(i1) = 0._RKIND
          end do
       end if
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if

    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep1d_d

  subroutine tr_allocatep1d_ll(array1d,begin_dim1,end_dim1,var_name,ocategory)
    integer*8, dimension(:)  , pointer    :: array1d
    integer*8                , intent(in) :: begin_dim1
    integer*8                , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    integer   :: err
    integer*8 :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep1d_ll

  subroutine tr_allocatep1d_ld(array1d,begin_dim1,end_dim1,var_name,ocategory,pzeroing)
    real(RKIND), dimension(:), pointer    :: array1d
    integer*8                , intent(in) :: begin_dim1
    integer*8                , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    logical        , optional, intent(in) :: pzeroing
    integer        , optional, intent(in) :: ocategory
    logical :: zeroing
    integer :: category

    integer   :: err
    integer*8 :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (.not. present(pzeroing)) then
      zeroing = .true.
    else
      zeroing = pzeroing
    end if
    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array1D',var_name,category)
    if (err.eq.0) then
       if (zeroing) then
          do i1 = begin_dim1,end_dim1
             array1d(i1) = 0._RKIND
          end do
       end if
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if

    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep1d_ld

  subroutine tr_allocatep1d_il(array1d,begin_dim1,end_dim1,var_name,ocategory)
    integer*8, dimension(:)  , pointer    :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep1d_il

  subroutine tr_allocatep1d_li(array1d,begin_dim1,end_dim1,var_name,ocategory)
    integer, dimension(:)    , pointer    :: array1d
    integer*8                , intent(in) :: begin_dim1
    integer*8                , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    integer   :: err
    integer*8 :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep1d_li

  !---------------------------------------- 
  ! memory allocation for a 2D array
  !----------------------------------------
  subroutine tr_allocatep2d_i(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name,ocategory)
    integer    , dimension(:,:), pointer    :: array2d
    integer                    , intent(in) :: begin_dim1
    integer                    , intent(in) :: end_dim1
    integer                    , intent(in) :: begin_dim2
    integer                    , intent(in) :: end_dim2
    character*(*)  , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array2D',var_name,category)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = 0
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep2d_i

  subroutine tr_allocatep2d_d(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name,ocategory,pzeroing)
    real(RKIND)  , dimension(:,:), pointer    :: array2d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1  
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    character*(*)    , intent(in)             :: var_name
    integer            , optional, intent(in) :: ocategory
    logical            , optional, intent(in) :: pzeroing
    logical :: zeroing
    integer :: category

    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (.not. present(pzeroing)) then
      zeroing = .true.
    else
      zeroing = pzeroing
    end if
    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array2D',var_name,category)
    if (err.eq.0) then
       if (zeroing) then
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array2d(i1,i2) = 0._RKIND
             end do
          end do
       end if
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep2d_d

  subroutine tr_allocatep2d_ll(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name,ocategory)
    integer*8  , dimension(:,:), pointer    :: array2d
    integer*8                  , intent(in) :: begin_dim1
    integer*8                  , intent(in) :: end_dim1
    integer*8                  , intent(in) :: begin_dim2
    integer*8                  , intent(in) :: end_dim2
    character*(*)  , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array2D',var_name,category)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = 0
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep2d_ll

  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine tr_allocatep3d_d(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name,ocategory,pzeroing)
    real(RKIND), dimension(:,:,:), pointer    :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name
    integer            , optional, intent(in) :: ocategory
    logical            , optional, intent(in) :: pzeroing
    logical :: zeroing
    integer category
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (.not. present(pzeroing)) then
      zeroing = .true.
    else
      zeroing = pzeroing
    end if

    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* sizeof(myreal)
    call tr_memwriteadd(size_array,'double array3D',var_name,category)
    if (err.eq.0) then
       if (zeroing) then
          do i3 = begin_dim3,end_dim3
             do i2 = begin_dim2,end_dim2
                do i1 = begin_dim1,end_dim1
                   array3d(i1,i2,i3) = 0._RKIND
                end do
             end do
          end do
       end if
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep3d_d

  !---------------------------------------- 
  ! memory allocation for a 4D array
  !----------------------------------------
  subroutine tr_allocatep4d_d(array4d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,begin_dim4,end_dim4,&
       var_name,ocategory,pzeroing)
    real(RKIND), dimension(:,:,:,:), pointer    :: array4d
    integer                        , intent(in) :: begin_dim1
    integer                        , intent(in) :: end_dim1
    integer                        , intent(in) :: begin_dim2
    integer                        , intent(in) :: end_dim2
    integer                        , intent(in) :: begin_dim3
    integer                        , intent(in) :: end_dim3
    integer                        , intent(in) :: begin_dim4
    integer                        , intent(in) :: end_dim4
    character*(*)                  , intent(in) :: var_name
    integer              , optional, intent(in) :: ocategory
    logical              , optional, intent(in) :: pzeroing
    logical   :: zeroing
    integer   :: category
    integer   :: err
    integer   :: i1, i2, i3, i4
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (.not. present(pzeroing)) then
      zeroing = .true.
    else
      zeroing = pzeroing
    end if

    allocate(array4d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3,begin_dim4:end_dim4),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)*(end_dim4-begin_dim4+1)* sizeof(myreal)
    call tr_memwriteadd(size_array,'double array4D',var_name,category)
    if (err.eq.0) then
       if (zeroing) then
          do i4 = begin_dim4,end_dim4
             do i3 = begin_dim3,end_dim3
                do i2 = begin_dim2,end_dim2
                   do i1 = begin_dim1,end_dim1
                      array4d(i1,i2,i3,i4) = 0._RKIND
                   end do
                end do
             end do
          end do
       end if
    else
       print *,'problem in allocating ',var_name
       print *,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep4d_d

  subroutine tr_allocatep3d_i(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name,ocategory)
    integer, dimension(:,:,:)    , pointer    :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* sizeof(myint)
    call tr_memwriteadd(size_array,'integer array3D',var_name,category)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = 0
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep3d_i

  subroutine tr_allocatep4d_i(array4d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3, begin_dim4, end_dim4, &
       var_name,ocategory)
    integer, dimension(:,:,:,:)  , pointer    :: array4d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    integer                      , intent(in) :: begin_dim4
    integer                      , intent(in) :: end_dim4
    character*(*)    , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category

    integer   :: err
    integer   :: i1, i2, i3, i4
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array4d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3,begin_dim4:end_dim4),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* &
         (end_dim4-begin_dim4+1)* sizeof(myint)
    call tr_memwriteadd(size_array,' integer array4D',var_name,category)
    if (err.eq.0) then
       do i4 = begin_dim4, end_dim4
          do i3 = begin_dim3,end_dim3
             do i2 = begin_dim2,end_dim2
                do i1 = begin_dim1,end_dim1
                   array4d(i1,i2,i3,i4) = 0
                end do
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocatep4d_i

  !---------------------------------------- 
  ! memory allocation for a 1D array
  !----------------------------------------
  subroutine tr_allocate1d_ch(array1d,begin_dim1,end_dim1,var_name,ocategory)
    character, dimension(:)    , allocatable :: array1d
    integer                    , intent(in) :: begin_dim1
    integer                    , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = ''
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate1d_ch

  subroutine tr_allocate1d_i(array1d,begin_dim1,end_dim1,var_name,ocategory)
    integer, dimension(:)    , allocatable :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory


    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate1d_i

  subroutine tr_allocate1d_d(array1d,begin_dim1,end_dim1,var_name,ocategory,pzeroing)
    real(RKIND), dimension(:), allocatable :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    logical, optional, intent(in) :: pzeroing
    logical :: zeroing
    integer category
    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (.not. present(pzeroing)) then
      zeroing = .true.
    else
      zeroing = pzeroing
    end if

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array1D',var_name,category)
    if (err.eq.0) then
       if (zeroing) then
          do i1 = begin_dim1,end_dim1
             array1d(i1) = 0._RKIND
          end do
       end if
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if

    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate1d_d

  subroutine tr_allocate1d_c(array1d,begin_dim1,end_dim1,var_name,ocategory)
    complex      , dimension(:), allocatable :: array1d
    integer                     , intent(in) :: begin_dim1
    integer                     , intent(in) :: end_dim1
    character*(*)   , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category

    integer   :: err
    integer   :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(mycomp)
    call tr_memwriteadd(size_array,'complex array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = cmplx(0,0)
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate1d_c

  subroutine tr_allocate1d_ll(array1d,begin_dim1,end_dim1,var_name,ocategory)
    integer*8, dimension(:)  , allocatable :: array1d
    integer*8                , intent(in) :: begin_dim1
    integer*8                , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory


    integer   :: err
    integer*8 :: i1
    integer*8 :: size_array 
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate1d_ll

  subroutine tr_allocate1d_ld(array1d,begin_dim1,end_dim1,var_name,ocategory,pzeroing)
    real(RKIND), dimension(:), allocatable :: array1d
    integer*8                , intent(in) :: begin_dim1
    integer*8                , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    logical, optional, intent(in) :: pzeroing
    logical :: zeroing
    integer category
    integer   :: err
    integer*8 :: i1
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (.not. present(pzeroing)) then
      zeroing = .true.
    else
      zeroing = pzeroing
    end if

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array1D',var_name,category)
    if (err.eq.0) then
       if (zeroing) then
          do i1 = begin_dim1,end_dim1
             array1d(i1) = 0._RKIND
          end do
       end if
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if

    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate1d_ld

  subroutine tr_allocate1d_il(array1d,begin_dim1,end_dim1,var_name,ocategory)
    integer*8, dimension(:)  , allocatable :: array1d
    integer                  , intent(in) :: begin_dim1
    integer                  , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory


    integer   :: err
    integer   :: i1
    integer*8 :: size_array 
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate1d_il

  subroutine tr_allocate1d_li(array1d,begin_dim1,end_dim1,var_name,ocategory)
    integer, dimension(:)    , allocatable :: array1d
    integer*8                , intent(in) :: begin_dim1
    integer*8                , intent(in) :: end_dim1
    character*(*), intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory


    integer   :: err
    integer*8 :: i1
    integer*8 :: size_array 
    integer category
    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array1d(begin_dim1:end_dim1),stat=err)
    size_array = (end_dim1-begin_dim1+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array1D',var_name,category)
    if (err.eq.0) then
       do i1 = begin_dim1,end_dim1
          array1d(i1) = 0
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate1d_li


  !---------------------------------------- 
  ! memory allocation for a 2D array
  !----------------------------------------
  subroutine tr_allocate2d_i(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name,ocategory)
    integer    , dimension(:,:), allocatable :: array2d
    integer                    , intent(in) :: begin_dim1
    integer                    , intent(in) :: end_dim1
    integer                    , intent(in) :: begin_dim2
    integer                    , intent(in) :: end_dim2
    character*(*)  , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category
    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array2D',var_name,category)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = 0
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate2d_i

  subroutine tr_allocate2d_d(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name,ocategory,pzeroing)
    real(RKIND)  , dimension(:,:), allocatable :: array2d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1  
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    character*(*)    , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory
    logical, optional, intent(in) :: pzeroing
    logical :: zeroing
    integer category
    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (.not. present(pzeroing)) then
      zeroing = .true.
    else
      zeroing = pzeroing
    end if

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myreal)
    call tr_memwriteadd(size_array,'double array2D',var_name,category)
    if (err.eq.0) then
       if (zeroing) then
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array2d(i1,i2) = 0._RKIND
             end do
          end do
       end if
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate2d_d

  subroutine tr_allocate2d_c(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name,ocategory)
    complex      , dimension(:,:), allocatable :: array2d
    integer                       , intent(in) :: begin_dim1
    integer                       , intent(in) :: end_dim1
    integer                       , intent(in) :: begin_dim2
    integer                       , intent(in) :: end_dim2
    character*(*)     , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category
    integer   :: err
    integer   :: i1, i2
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(mycomp)
    call tr_memwriteadd(size_array,'complex array2D',var_name,category)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = cmplx(0,0)
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate2d_c

  subroutine tr_allocate2d_ll(array2d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,var_name,ocategory)
    integer*8  , dimension(:,:), allocatable :: array2d
    integer*8                  , intent(in) :: begin_dim1
    integer*8                  , intent(in) :: end_dim1
    integer*8                  , intent(in) :: begin_dim2
    integer*8                  , intent(in) :: end_dim2
    character*(*)  , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category
    integer   :: err
    integer*8 :: i1, i2
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array2d(begin_dim1:end_dim1,begin_dim2:end_dim2), &
         stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1) * sizeof(myint)
    call tr_memwriteadd(size_array,'integer array2D',var_name,category)
    if (err.eq.0) then
       do i2 = begin_dim2,end_dim2
          do i1 = begin_dim1,end_dim1
             array2d(i1,i2) = 0
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate2d_ll


  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine tr_allocate3d_i(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name,ocategory)
    integer    , dimension(:,:,:), allocatable :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* SIZEOF(myint)
    call tr_memwriteadd(size_array,'integer array3D',var_name,category)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = 0
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate3d_i

  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine tr_allocate3d_c(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name,ocategory)
    complex   , dimension(:,:,:), allocatable :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* SIZEOF(mycomp)
    call tr_memwriteadd(size_array,'integer array3D',var_name,category)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = cmplx(0,0)
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate3d_c

  !---------------------------------------- 
  ! memory allocation for a 3D array
  !----------------------------------------
  subroutine tr_allocate3d_d(array3d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,var_name,ocategory)
    real(RKIND), dimension(:,:,:), allocatable :: array3d
    integer                      , intent(in) :: begin_dim1
    integer                      , intent(in) :: end_dim1
    integer                      , intent(in) :: begin_dim2
    integer                      , intent(in) :: end_dim2
    integer                      , intent(in) :: begin_dim3
    integer                      , intent(in) :: end_dim3
    character*(*)    , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category
    integer   :: err
    integer   :: i1, i2, i3
    integer*8 :: size_array 


    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array3d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* sizeof(myreal)
    call tr_memwriteadd(size_array,'double array3D',var_name,category)
    if (err.eq.0) then
       do i3 = begin_dim3,end_dim3
          do i2 = begin_dim2,end_dim2
             do i1 = begin_dim1,end_dim1
                array3d(i1,i2,i3) = 0._RKIND
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate3d_d

  !---------------------------------------- 
  ! memory allocation for a 4D array
  !----------------------------------------
  subroutine tr_allocate4d_d(array4d,begin_dim1,end_dim1, &
       begin_dim2,end_dim2,begin_dim3,end_dim3,begin_dim4,end_dim4, &
       var_name,ocategory)
    real(RKIND), dimension(:,:,:,:), allocatable :: array4d
    integer                        , intent(in) :: begin_dim1
    integer                        , intent(in) :: end_dim1
    integer                        , intent(in) :: begin_dim2
    integer                        , intent(in) :: end_dim2
    integer                        , intent(in) :: begin_dim3
    integer                        , intent(in) :: end_dim3
    integer                        , intent(in) :: begin_dim4
    integer                        , intent(in) :: end_dim4
    character*(*)    , intent(in)             :: var_name
    integer  , optional, intent(in) :: ocategory

    integer category
    integer   :: err
    integer   :: i1, i2, i3, i4
    integer*8 :: size_array 

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    allocate(array4d(begin_dim1:end_dim1,begin_dim2:end_dim2, &
         begin_dim3:end_dim3,begin_dim4:end_dim4),stat=err)
    size_array = (end_dim1-begin_dim1+1) * &
         (end_dim2-begin_dim2+1)* &
         (end_dim3-begin_dim3+1)* &
         (end_dim4-begin_dim4+1)* sizeof(myreal)
    call tr_memwriteadd(size_array,'double array4D',var_name,category)
    if (err.eq.0) then
       do i4 = begin_dim4,end_dim4
          do i3 = begin_dim3,end_dim3
             do i2 = begin_dim2,end_dim2
                do i1 = begin_dim1,end_dim1
                   array4d(i1,i2,i3,i4) = 0._RKIND
                end do
             end do
          end do
       end do
    else
       print*,'problem in allocating ',var_name
       print*,'-> required memory (in Bytes) = ',size_array
       stop
    end if
    nb_allocate(category) = nb_allocate(category) + size_array
    max_allocate = max(max_allocate,SUM(nb_allocate(MIN_CAT:MAX_CAT)))
  end subroutine tr_allocate4d_d

  !---------------------------------------- 
  ! memory deallocation of array 1D
  !----------------------------------------
  subroutine tr_deallocatep1d_i(array1d,var_name,ocategory)
    integer, dimension(:), pointer     :: array1d
    character*(*)        , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (associated(array1d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array1d)
       call tr_memwritedel(var_name,category)
       deallocate(array1d)
       array1d => NULL()
    end if
  end subroutine tr_deallocatep1d_i

  subroutine tr_deallocatep1d_d(array1d,var_name,ocategory)
    real(RKIND), dimension(:), pointer :: array1d
    character*(*)        , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (associated(array1d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array1d)
       call tr_memwritedel(var_name,category)
       deallocate(array1d)
       array1d => NULL()
    end if
  end subroutine tr_deallocatep1d_d

  subroutine tr_deallocatep1d_l(array1d,var_name,ocategory)
    integer*8, dimension(:), pointer     :: array1d
    character*(*)        , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if
    if (associated(array1d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array1d)
       call tr_memwritedel(var_name,category)
       deallocate(array1d)
       array1d => NULL()
    end if
  end subroutine tr_deallocatep1d_l

  !---------------------------------------- 
  ! memory deallocation of array 2D
  !----------------------------------------
  subroutine tr_deallocatep2d_i(array2d,var_name,ocategory)
    integer, dimension(:,:) , pointer  :: array2d
    character*(*)           , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (associated(array2d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array2d)
       call tr_memwritedel(var_name,category)
       deallocate(array2d)
       array2d => null()
    end if
  end subroutine tr_deallocatep2d_i

  subroutine tr_deallocatep2d_d(array2d,var_name,ocategory)
    real(RKIND), dimension(:,:), pointer :: array2d
    character*(*)              , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (associated(array2d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array2d)
       call tr_memwritedel(var_name,category)
       deallocate(array2d)
       array2d => null()
    end if
  end subroutine tr_deallocatep2d_d

  subroutine tr_deallocatep2d_l(array2d,var_name,ocategory)
    integer*8, dimension(:,:) , pointer  :: array2d
    character*(*)           , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (associated(array2d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array2d)
       call tr_memwritedel(var_name,category)
       deallocate(array2d)
       array2d => null()
    end if
  end subroutine tr_deallocatep2d_l

  !---------------------------------------- 
  ! memory deallocation of array 3D
  !----------------------------------------
  subroutine tr_deallocatep3d_d(array3d,var_name,ocategory)
    real(RKIND), dimension(:,:,:) , pointer :: array3d
    character*(*)                 , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (associated(array3d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array3d)
       call tr_memwritedel(var_name,category)
       deallocate(array3d)
       array3d => null()
    end if
  end subroutine tr_deallocatep3d_d

  subroutine tr_deallocatep3d_i(array3d,var_name,ocategory)
    integer, dimension(:,:,:) , pointer :: array3d
    character*(*)                 , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (associated(array3d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array3d)
       call tr_memwritedel(var_name,category)
       deallocate(array3d)
       array3d => null()
    end if
  end subroutine tr_deallocatep3d_i

  !---------------------------------------- 
  ! memory deallocation of array 4D
  !----------------------------------------
  subroutine tr_deallocatep4d_d(array4d,var_name,ocategory)
    real(RKIND), dimension(:,:,:,:) , pointer :: array4d
    character*(*)                   , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (associated(array4d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array4d)
       call tr_memwritedel(var_name,category)
       deallocate(array4d)
       array4d => null()
    end if
  end subroutine tr_deallocatep4d_d

  subroutine tr_deallocatep4d_i(array4d,var_name,ocategory)
    integer, dimension(:,:,:,:)   , pointer :: array4d
    character*(*)                 , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (associated(array4d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array4d)
       call tr_memwritedel(var_name,category)
       deallocate(array4d)
       array4d => null()
    end if
  end subroutine tr_deallocatep4d_i


  !---------------------------------------- 
  ! memory deallocation of array 1D
  !----------------------------------------
  subroutine tr_deallocate1d_ch(array1d,var_name,ocategory)
    character, dimension(:), allocatable  :: array1d
    character*(*)          , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array1d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array1d)
       call tr_memwritedel(var_name,category)
       deallocate(array1d)
    end if
  end subroutine tr_deallocate1d_ch

  subroutine tr_deallocate1d_i(array1d,var_name,ocategory)
    integer, dimension(:), allocatable  :: array1d
    character*(*)        , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array1d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array1d)
       call tr_memwritedel(var_name,category)
       deallocate(array1d)
    end if
  end subroutine tr_deallocate1d_i

  subroutine tr_deallocate1d_d(array1d,var_name,ocategory)
    real(RKIND), dimension(:), allocatable :: array1d
    character*(*)        , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array1d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array1d)
       call tr_memwritedel(var_name,category)
       deallocate(array1d)
    end if
  end subroutine tr_deallocate1d_d

  subroutine tr_deallocate1d_c(array1d,var_name,ocategory)
    complex      , dimension(:),  allocatable :: array1d
    character*(*)               , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array1d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array1d)
       call tr_memwritedel(var_name,category)
       deallocate(array1d)
    end if
  end subroutine tr_deallocate1d_c

  subroutine tr_deallocate1d_l(array1d,var_name,ocategory)
    integer*8, dimension(:), allocatable  :: array1d
    character*(*)        , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array1d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array1d)
       call tr_memwritedel(var_name,category)
       deallocate(array1d)
    end if
  end subroutine tr_deallocate1d_l

  !---------------------------------------- 
  ! memory deallocation of array 2D
  !----------------------------------------
  subroutine tr_deallocate2d_i(array2d,var_name,ocategory)
    integer, dimension(:,:) , allocatable  :: array2d
    character*(*)           , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array2d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array2d)
       call tr_memwritedel(var_name,category)
       deallocate(array2d)
    end if
  end subroutine tr_deallocate2d_i

  subroutine tr_deallocate2d_d(array2d,var_name,ocategory)
    real(RKIND), dimension(:,:), allocatable :: array2d
    character*(*)              , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array2d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array2d)
       call tr_memwritedel(var_name,category)
       deallocate(array2d)
    end if
  end subroutine tr_deallocate2d_d

  subroutine tr_deallocate2d_c(array2d,var_name,ocategory)
    complex      , dimension(:,:), allocatable :: array2d
    character*(*)                 , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array2d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array2d)
       call tr_memwritedel(var_name,category)
       deallocate(array2d)
    end if
  end subroutine tr_deallocate2d_c

  subroutine tr_deallocate2d_l(array2d,var_name,ocategory)
    integer*8, dimension(:,:) , allocatable  :: array2d
    character*(*)           , intent(in)  :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array2d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array2d)
       call tr_memwritedel(var_name,category)
       deallocate(array2d)
    end if
  end subroutine tr_deallocate2d_l

  !---------------------------------------- 
  ! memory deallocation of array 3D
  !----------------------------------------
  subroutine tr_deallocate3d_d(array3d,var_name,ocategory)
    real(RKIND), dimension(:,:,:) , allocatable :: array3d
    character*(*)                 , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array3d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array3d)
       call tr_memwritedel(var_name,category)
       deallocate(array3d)
    end if
  end subroutine tr_deallocate3d_d

  subroutine tr_deallocate4d_d(array4d,var_name,ocategory)
    real(RKIND), dimension(:,:,:,:) , allocatable :: array4d
    character*(*)                   , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array4d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array4d)
       call tr_memwritedel(var_name,category)
       deallocate(array4d)
    end if
  end subroutine tr_deallocate4d_d

  subroutine tr_deallocate3d_c(array3d,var_name,ocategory)
    complex      , dimension(:,:,:), allocatable :: array3d
    character*(*)                   , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array3d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array3d)
       call tr_memwritedel(var_name,category)
       deallocate(array3d)
    end if
  end subroutine tr_deallocate3d_c

  subroutine tr_deallocate3d_i(array3d,var_name,ocategory)
    integer      , dimension(:,:,:), allocatable :: array3d
    character*(*)                   , intent(in) :: var_name
    integer  , optional, intent(in) :: ocategory
    integer category

    if (.not. present(ocategory)) then
      category = CAT_UNKNOWN
    else
      category = ocategory
    end if

    if (allocated(array3d)) then
       nb_allocate(category) = nb_allocate(category) - sizeof(array3d)
       call tr_memwritedel(var_name,category)
       deallocate(array3d)
    end if
  end subroutine tr_deallocate3d_i


  !***********************************************
  !  function for program analysis
  !***********************************************
  subroutine tr_print_memsize(label)
    implicit none
    character(*), intent(in) :: label
    integer*8, parameter :: GBconst = 1024_8*1024_8*1024_8
    integer*8, parameter :: MBconst = 1024_8*1024_8
    integer*8, parameter :: KBconst = 1024_8
    integer :: uout, j
    integer*8 :: scount, dcount, rcount, lcount, pcount
    
    rcount = KBconst * get_memory_inkb("VmRSS")
    open(uout_mem, file = trace_file, status = 'OLD', &
         position = 'APPEND', form = 'FORMATTED')
#ifdef MEMTRACE
    do j = MIN_CAT,MAX_CAT
      write(uout_mem,'(A20,A37,A8,A5,1f10.3,A)') label, &
        'memsize allocated within Jorek (cat=',cat_name(j),') = ', &
        nb_allocate(j)/dfloat(MBconst), ' MBytes'
    end do
#endif
    if (SUM(nb_allocate(MIN_CAT:MAX_CAT)).gt.GBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)') label,&
            'memsize allocated within Jorek (total) = ', &
            SUM(nb_allocate(MIN_CAT:MAX_CAT))/dfloat(GBconst), ' GBytes'
    else 
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize allocated within Jorek (total) = ', &
            SUM(nb_allocate(MIN_CAT:MAX_CAT))/dfloat(MBconst), ' MBytes'
    end if
    lcount = rcount - SUM(nb_allocate(MIN_CAT:MAX_CAT))
    if (lcount .lt. 0) then
       write(uout_mem,'(A20,A50,1G10.3,A)') label, &
            'memsize occupied by libraries/others = ', &
            dfloat(lcount), ' Bytes'
    else if (lcount.gt.GBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize occupied by libraries/others = ', &
            lcount/dfloat(GBconst), ' GBytes'
    else if (lcount.gt.MBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize occupied by libraries/others = ', &
            lcount/dfloat(MBconst), ' MBytes'
    else if (lcount .gt. KBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize occupied by libraries/others = ', &
            lcount/dfloat(KBconst), ' KBytes'
    else
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize occupied by libraries/others = ', &
            lcount, ' Bytes'
    end if
#ifdef USE_PASTIX
    call pastix_getmem(pcount)
    if (pcount .gt. GBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize PaStiX  = ', &
            pcount/dfloat(GBconst), ' GBytes'
    else if (pcount.gt.MBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize PaStiX = ', &
            pcount/dfloat(MBconst), ' MBytes'
    else
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize PaStiX = ', &
            pcount/dfloat(KBconst), ' KBytes'
    end if
#endif    
    if (rcount.gt.GBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize total   (RSS) = ', &
            rcount/dfloat(GBconst), ' GBytes'
    else if (rcount.gt.MBconst) then
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize total   (RSS) = ', &
            rcount/dfloat(MBconst), ' MBytes'
    else
       write(uout_mem,'(A20,A50,1f10.3,A)') label, &
            'memsize total   (RSS) = ', &
            rcount/dfloat(KBconst), ' KBytes'
    end if
    close(uout_mem)
  end subroutine tr_print_memsize
end module tr_module
