!> Contains data structures and some routines related to the global matrix and right hand sides.
!!
!! The sparse matrix and the corresponding right hand side represent the system of equations
!! given to the solver for calculating the time evolution of the physical quantities.
module global_distributed_matrix
  
  use mod_integer_types

  implicit none
  
  public
  
  ! --- The complex harmonic matrix 
#ifdef USE_COMPLEX_PRECOND
  double complex,        allocatable, target :: A_cmplx(:)          !< Distributed harmonic matrix
  double complex,        allocatable, target :: rhs_cmplx(:)        !< Distributed harmonic right hand side
  double complex,        allocatable, target :: rhs_cmplx_guess(:)  !< Guess solution for GMRES
  double complex,        allocatable, target :: rhs_cmplx_sol(:)    !< Solution from GMRES
  integer(kind=int_all), allocatable, target :: irn_cmplx(:)        !< Row indices for coordinate format sparse matrix (or CSR)
  integer(kind=int_all), allocatable, target :: jcn_cmplx(:)        !< Column indices for coordinate format sparse matrix (or CSR)
  integer(kind=int_all)                      :: n_cmplx, nz_cmplx                       
#endif
 
  contains
  
  
  
  !> Determine the matrix row or column for given values of ::i_index, ::i_var, and ::i_tor.
  integer pure function det_row_col(i_index, i_var, i_tor, i_tor_min, i_tor_max)
    
    use mod_parameters, only: n_tor, n_var
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_index   !< node%index property for the node and node degree of freedom
    integer, intent(in) :: i_var     !< Variable number
    integer, intent(in) :: i_tor     !< Toroidal mode number
    integer, intent(in) :: i_tor_min !< Minimum Toroidal mode number
    integer, intent(in) :: i_tor_max !< Maximum Toroidal mode number

    ! --- Local variables
    integer :: n_tor_local           ! local toroidal number

    n_tor_local = i_tor_max - i_tor_min + 1 
    det_row_col =  n_tor_local * n_var * (i_index-1) + n_tor_local * (i_var-1) + i_tor - i_tor_min + 1
    
  end function det_row_col
  
  
  
  
  
  
  !> Determine the position of a matrix entry given by its row and column positions (::i_row and
  !! ::j_col) in the sparse matrix structure.
  integer(kind=int_all) pure function det_sparse_pos(i_row, j_col, index_min, a_mat)
    
    use mod_integer_types
    use data_structure, only: type_SP_MATRIX

    implicit none
    
    ! --- Routine parameters
    integer,                            intent(in)    :: i_row                   !< Matrix row
    integer,                            intent(in)    :: j_col                   !< Matrix column
    integer,                            intent(in)    :: index_min               !< Smallest block index dealt with by current MPI proc
    type(type_SP_MATRIX), intent(in)                  :: a_mat
    
    ! --- Local variables
    integer :: i_block, j_block           ! Block indices
    integer :: i_row_block, j_col_block   ! Row and column in the block
    integer :: ij_sparse_block            ! Position of the block in the sparse matrix structure
    integer :: i_block_local              ! Block index at local MPI proc
    integer(kind=int_all) :: i
    
    i_block       = (i_row-1) / a_mat%block_size + 1
    i_block_local = i_block - index_min + 1
    i_row_block   = i_row - (i_block-1) * a_mat%block_size
    
    j_block       = (j_col-1) / a_mat%block_size + 1
    j_col_block   = j_col - (j_block-1) * a_mat%block_size
    
    ! --- Determine the position of the block in the sparse matrix.
    det_sparse_pos = -999999999 ! return a negative value if not found
    do i = 1, a_mat%ijA_size(i_block_local)
      if (a_mat%irn_jcn(i_block_local,i) == j_block ) then ! Block index found?
        ij_sparse_block = a_mat%ijA_index(i_block_local,i)
        det_sparse_pos  = ij_sparse_block-1 + a_mat%block_size*(i_row_block-1) + j_col_block
        exit
      end if
    end do
    
  end function det_sparse_pos
  
  
  
  
  
  
  !!> Determine the position of a matrix block given by its row and column positions
  !!! in the sparse matrix structure.
  !integer(kind=int_all) pure function det_sparse_pos_block(i_row, j_col, index_min, a_mat)
  !  use data_structure, only: type_SP_MATRIX
  !  
  !  implicit none
  !  
  !  ! --- Routine parameters
  !  integer, intent(in) :: i_row     !< Matrix row
  !  integer, intent(in) :: j_col     !< Matrix column
  !  integer, intent(in) :: index_min !< Smallest block index dealt with by current MPI proc
  !  type(type_SP_MATRIX), intent(in) :: a_mat
  !  
  !  ! --- Local variables
  !  integer :: i_block, j_block           ! Block indices
  !  integer :: i_block_local              ! Block index at local MPI proc
  !  integer :: i
  !  
  !  i_block       = (i_row-1) / a_mat%block_size + 1
  !  i_block_local = i_block - index_min + 1
  !  j_block       = (j_col-1) / a_mat%block_size + 1
  !  
  !  ! --- Determine the position of the block in the sparse matrix.
  !  det_sparse_pos_block = -999999999 !###
  !  do i = 1, ijA_size(i_block_local)
  !    if ( irn_jcn(i_block_local,i) == j_block ) then ! Block index found?
  !      det_sparse_pos_block = ijA_index(i_block_local,i)
  !      exit
  !    end if
  !  end do
  !  
  !end function det_sparse_pos_block
  
  
  
  
  
  
  !> Initialize the (freeboundary related) row and column numbers in the sparse matrix structure.
  subroutine global_matrix_structure_vacuum(node_list, bnd_node_list, a_mat, i_tor_min, i_tor_max)
    
    use mod_parameters, only: n_tor, n_var
    use data_structure, only: type_node_list, type_bnd_node_list, type_SP_MATRIX
    use mod_integer_types
    
    implicit none
    
    ! --- Routine parameters
    type(type_node_list),               intent(in)    :: node_list            !< List of grid nodes
    type(type_bnd_node_list),           intent(in)    :: bnd_node_list        !< List of boundary grid nodes
    integer,                            intent(in)    :: i_tor_min, i_tor_max !< Toroidal mode numbers 
    type(type_SP_MATRIX)                              :: a_mat
    ! --- Local variables
    integer :: l_node_bnd, l_dof, l_node, l_dir, l_index, l_tor, l_var, l_row
    integer :: j_node_bnd, j_dof, j_node, j_dir, j_index, j_tor, j_var, j_col
    integer :: sparsepos
    
    !$omp parallel do                                                                    &
    !$omp default(none)                                                                  &
    !$omp shared(bnd_node_list, node_list, a_mat, i_tor_min, i_tor_max)                  &
    !$omp private(l_node_bnd, l_dof, l_tor, l_var, j_node_bnd, j_dof,  j_tor,            &
    !$omp         j_var, l_node, l_dir, l_index, l_row, j_node, j_dir, j_index,          &
    !$omp         j_col, sparsepos)                                                      &
    !$omp schedule(dynamic,1)
    ! --- Select a matrix row
    do l_node_bnd = 1, bnd_node_list%n_bnd_nodes
      do l_dof = 1, 2
        l_node      = bnd_node_list%bnd_node(l_node_bnd)%index_jorek
        l_dir       = bnd_node_list%bnd_node(l_node_bnd)%direction(l_dof)
        l_index     = node_list%node(l_node)%index(l_dir)
        if ( (l_index < a_mat%my_ind_min) .or. (l_index > a_mat%my_ind_max) ) cycle ! Is the current MPI thread in charge?
        do l_tor = i_tor_min, i_tor_max 
          do l_var = 1, n_var
            l_row = det_row_col(l_index, l_var, l_tor, i_tor_min, i_tor_max)
            
            ! --- Select a matrix column
            do j_node_bnd = 1, bnd_node_list%n_bnd_nodes
              do j_dof = 1, 2
                j_node      = bnd_node_list%bnd_node(j_node_bnd)%index_jorek
                j_dir       = bnd_node_list%bnd_node(j_node_bnd)%direction(j_dof)
                j_index     = node_list%node(j_node)%index(j_dir)
                do j_tor = i_tor_min, i_tor_max
                  do j_var = 1, n_var
                    j_col = det_row_col(j_index, j_var, j_tor, i_tor_min, i_tor_max)
                    
                    ! --- Determine which position in the sparse matrix data structure corresponds
                    !     to the matrix entry at l_row, j_col.
                    sparsepos = det_sparse_pos(l_row, j_col, a_mat%my_ind_min, a_mat)
                    
                    ! --- Set row and column numbers in the sparse matrix data structure
                    a_mat%irn(sparsepos) = l_row
                    a_mat%jcn(sparsepos) = j_col
                    
                  end do
                end do
              end do
            end do
            
          end do
        end do
      end do
    end do
    !$omp end parallel do
    
  end subroutine global_matrix_structure_vacuum
  
  
  
end module global_distributed_matrix
