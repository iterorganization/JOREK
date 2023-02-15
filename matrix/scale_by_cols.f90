!> Performs column scaling of sparse matrix  
  subroutine scale_by_cols(a_mat)
    use data_structure, only: type_SP_MATRIX, type_RHS
    use mpi_mod
    use mod_integer_types
    
    implicit none
    
    type(type_SP_MATRIX)    :: a_mat
    type(type_RHS)          :: lcs     ! local column scaling
    integer(kind=int_all)   :: j, k
    integer                 :: counts, ierr
    
    if (associated(a_mat%column_scaling)) then
      deallocate(a_mat%column_scaling); a_mat%column_scaling => Null()
    endif
    allocate(a_mat%column_scaling(a_mat%ng))
    lcs%n = a_mat%ng
    allocate(lcs%val(lcs%n))

    a_mat%column_scaling(1:a_mat%ng) = 1.d-20
    lcs%val(1:lcs%n) = 1.d-20
    do k=1,a_mat%nnz
      j = a_mat%jcn(k)
      lcs%val(j) = max(lcs%val(j),abs(a_mat%val(k)))
    enddo

    counts = a_mat%ng
    call MPI_AllReduce(lcs%val,a_mat%column_scaling,counts,MPI_DOUBLE_PRECISION,MPI_MAX,a_mat%comm,ierr)
    
    do k = 1, a_mat%nnz
      j = a_mat%jcn(k)
      a_mat%val(k) = a_mat%val(k)/a_mat%column_scaling(j)
    enddo
    a_mat%scaled = .true.
    
    deallocate(lcs%val)
    return
    
  end subroutine scale_by_cols