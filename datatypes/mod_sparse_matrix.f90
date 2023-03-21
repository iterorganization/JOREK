module mod_sparse_matrix
  use mod_integer_types
  implicit none
  
  private
  public :: type_SP_MATRIX

  !> Sparse matrix type (generally distributed)
  type type_SP_MATRIX
    integer(kind=int_all), dimension(:), pointer   :: irn => Null()
    integer(kind=int_all), dimension(:), pointer   :: jcn => Null()
    real(kind=8), dimension(:), pointer            :: val => Null()
    integer(kind=int_all), dimension(:), pointer   :: ijA_size => Null()
    integer(kind=int_all), dimension(:,:), pointer :: ijA_index => Null()
    integer(kind=int_all), dimension(:,:), pointer :: irn_jcn => Null()
    
    real(kind=8), dimension(:), pointer          :: column_scaling => Null()    !< global column scaling, vector size of ng
    integer                                      :: indexing = 1         !< matrix indexing (1 is standart FORTRAN)
    integer(kind=int_all)                        :: ng = 0               !< matrix total rank
    integer(kind=int_all)                        :: nr = 0               !< number of local rows
    integer(kind=int_all)                        :: nc = 0               !< number of local cols
    integer(kind=int_all)                        :: nnz = 0              !< number of local nonzero entries
    integer, dimension(:), pointer               :: index_min => Null()  !< minimum node index in global range for all MPI ranks
    integer, dimension(:), pointer               :: index_max => Null()  !< maximum node index in global range for all MPI ranks
    integer                                      :: my_ind_min = 0
    integer                                      :: my_ind_max = 0
    integer                                      :: my_ind_size = 0
    integer(kind=int_all)                        :: maxsize = 0          !< maxval(ijA_size(:))
    integer                                      :: i_tor_min = 0        ! minimum toroidal Fourier number used in construction
    integer                                      :: i_tor_max = 0        ! maximum toroidal Fourier number used in construction
    integer                                      :: block_size = 1
    integer                                      :: comm = 0             !< communicator over which the matrix is distributed
    integer                                      :: ncpu = 1
    logical                                      :: scaled = .false.
    logical                                      :: row_distributed = .false.
    logical                                      :: col_distributed = .false.
    logical                                      :: reduced = .false. !< matrix is available on all comm ranks (not distribued)

  contains
    procedure :: copy_to
    procedure :: move_to
    procedure :: reset
  end type type_SP_MATRIX
  
  contains
  

! move one matrix structure into another
  subroutine move_to(self, mat_a, with_data)
    class(type_SP_MATRIX), intent(inout)    :: self
    class(type_SP_MATRIX), intent(inout)    :: mat_a
    logical                                 :: with_data

    if (with_data) then
      mat_a%irn            => self%irn; self%irn => null()
      mat_a%jcn            => self%jcn; self%jcn => null()
      mat_a%val            => self%val; self%val => null()

      mat_a%ijA_size       => self%ijA_size
      mat_a%ijA_index      => self%ijA_index
      mat_a%irn_jcn        => self%irn_jcn

      mat_a%column_scaling => self%column_scaling
      mat_a%index_min      => self%index_min
      mat_a%index_max      => self%index_max
    endif

    mat_a%indexing        = self%indexing
    mat_a%ng              = self%ng
    mat_a%nr              = self%nr
    mat_a%nc              = self%nc
    mat_a%nnz             = self%nnz
    mat_a%my_ind_min      = self%my_ind_min
    mat_a%my_ind_max      = self%my_ind_max
    mat_a%my_ind_size     = self%my_ind_size
    mat_a%maxsize         = self%maxsize
    mat_a%i_tor_min       = self%i_tor_min
    mat_a%i_tor_max       = self%i_tor_max
    mat_a%block_size      = self%block_size
    mat_a%comm            = self%comm
    mat_a%ncpu            = self%ncpu
    mat_a%scaled          = self%scaled
    mat_a%row_distributed = self%row_distributed
    mat_a%col_distributed = self%col_distributed
    mat_a%reduced         = self%reduced

    return
  end subroutine move_to

! copy one matrix structure into another
  subroutine copy_to(self, mat_a)
    class(type_SP_MATRIX), intent(inout)    :: self
    class(type_SP_MATRIX), intent(inout)    :: mat_a

    call mat_a%reset()

    call self%move_to(mat_a, with_data=.false.)

    if (associated(self%irn)) then
      allocate(mat_a%irn(mat_a%nnz))
      mat_a%irn(1:mat_a%nnz) = self%irn(1:self%nnz)
    endif
    if (associated(self%jcn)) then
      allocate(mat_a%jcn(mat_a%nnz))
      mat_a%jcn(1:mat_a%nnz) = self%jcn(1:self%nnz)
    endif
    if (associated(self%val)) then
      allocate(mat_a%val(mat_a%nnz))
      mat_a%val(1:mat_a%nnz) = self%val(1:self%nnz)
    endif
    if (associated(self%column_scaling)) then
      allocate(mat_a%column_scaling(mat_a%ng))
      mat_a%column_scaling(1:mat_a%ng) = self%column_scaling(1:self%ng)
    endif
    if (associated(self%index_min)) then
      allocate(mat_a%index_min(mat_a%ncpu))
      mat_a%index_min(1:mat_a%ncpu) = self%index_min(1:self%ncpu)
    endif
    if (associated(self%index_max)) then
      allocate(mat_a%index_max(mat_a%ncpu))
      mat_a%index_max(1:mat_a%ncpu) = self%index_max(1:self%ncpu)
    endif
    if (associated(self%ijA_size)) then
      allocate(mat_a%ijA_size(mat_a%my_ind_size))
      mat_a%ijA_size(1:mat_a%my_ind_size) = self%ijA_size(1:self%my_ind_size)
    endif
    if (associated(self%irn_jcn)) then
      allocate(mat_a%irn_jcn(mat_a%my_ind_size,mat_a%maxsize))
      mat_a%irn_jcn(1:mat_a%my_ind_size,1:mat_a%maxsize) = self%irn_jcn(1:self%my_ind_size,1:self%maxsize)
    endif
    if (associated(self%ijA_index)) then
      allocate(mat_a%ijA_index(mat_a%my_ind_size,mat_a%maxsize))
      mat_a%ijA_index(1:mat_a%my_ind_size,1:mat_a%maxsize) = self%ijA_index(1:self%my_ind_size,1:self%maxsize)
    endif    

    return
  end subroutine copy_to

  subroutine reset(self)
    class(type_SP_MATRIX), intent(inout)    :: self

    if (associated(self%irn)) then
      deallocate(self%irn); self%irn => Null()
    endif
    if (associated(self%jcn)) then
      deallocate(self%jcn); self%jcn => Null()
    endif
    if (associated(self%val)) then
      deallocate(self%val); self%val => Null()
    endif
    if (associated(self%ijA_size)) then
      deallocate(self%ijA_size); self%ijA_size => Null()
    endif
    if (associated(self%ijA_index)) then
      deallocate(self%ijA_index); self%ijA_index => Null()
    endif
    if (associated(self%irn_jcn)) then
      deallocate(self%irn_jcn); self%irn_jcn => Null()
    endif
    if (associated(self%column_scaling)) then
      deallocate(self%column_scaling); self%column_scaling => Null()
    endif
    if (associated(self%index_min)) then
      deallocate(self%index_min); self%index_min => Null()
    endif
    if (associated(self%index_max)) then
      deallocate(self%index_max); self%index_max => Null()
    endif

    self%indexing = 1
    self%ng = 0
    self%nr = 0
    self%nc = 0
    self%nnz = 0
    self%my_ind_min = 0
    self%my_ind_max = 0
    self%my_ind_size = 0
    self%maxsize = 0
    self%i_tor_min = 0
    self%i_tor_max = 0
    self%block_size = 1
    self%comm = 0
    self%scaled = .false.
    self%row_distributed = .false.
    self%col_distributed = .false.
    self%reduced = .false.
  end subroutine reset
  
end module mod_sparse_matrix