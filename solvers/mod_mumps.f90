module mod_mumps
#ifdef USE_MUMPS
  include "dmumps_struc.h"
  
  type type_MUMPS_SOLVER

    integer                                      :: comm = 0
    real(kind=8), dimension(:), pointer          :: solution_scaling => Null()    !< matrix column scaling to be applied to solution vector
    logical                                      :: initialized = .false.
    logical                                      :: analyzed    = .false.
    logical                                      :: equilibrium = .false.    
    logical                                      :: scaled      = .false.
    logical                                      :: refine      = .false.
    
    type(DMUMPS_STRUC)  :: mumps_par
    integer             :: mumps_ordering
    logical             :: use_BLR_compression
    real(kind=8)        :: epsilon_BLR
  
  end type type_MUMPS_SOLVER
  
  private
  public :: type_MUMPS_SOLVER, mumps_initialize, mumps_analyze, mumps_factorize, mumps_solve, mumps_finalize
  
  contains
  
  subroutine mumps_initialize(mmss,comm)
    use phys_module, only: mumps_ordering, epsilon_BLR, use_BLR_compression
    implicit none
    
    type(type_MUMPS_SOLVER)            :: mmss
    integer :: comm
    
    mmss%comm = comm ! normal communicator, doesnt work wor MPI_COMM_SELF
    
    mmss%mumps_par%COMM = mmss%comm                   ! Define a communicator for mumps
  
    mmss%mumps_par%JOB = -1
    mmss%mumps_par%SYM = 0
    mmss%mumps_par%PAR = 1
    mmss%mumps_par%icntl(13) = -1
  
    call  DMUMPS(mmss%mumps_par)                        ! Initialize an instance of mumps
    
    mmss%mumps_par%icntl(2)  = -1                       ! Verbosity levels
    mmss%mumps_par%icntl(3)  = -1
    mmss%mumps_par%icntl(4)  = 6
  
    mmss%mumps_par%icntl(14) = 20                       ! memory working space increase
    !mmss%mumps_par%icntl(15) = 0                       ! memory balance only, 1: flops      
  
    mmss%initialized = .true.
    
  end subroutine mumps_initialize
  
  subroutine mumps_analyze(mmss,a_mat)
  
    use data_structure, only: type_SP_MATRIX
    implicit none
    
    type(type_MUMPS_SOLVER) :: mmss
    type(type_SP_MATRIX)    :: a_mat
    
    
    mmss%mumps_par%A_loc   => a_mat%val
    mmss%mumps_par%irn_loc => a_mat%irn
    mmss%mumps_par%jcn_loc => a_mat%jcn
    
    mmss%mumps_par%n      = a_mat%ng
    mmss%mumps_par%nz_loc = a_mat%nnz    
    
    mmss%mumps_par%JOB = 1                                  ! Analysis, only needed when grid has changed
    
    mmss%mumps_par%icntl(7)  = mmss%mumps_ordering               ! ordering option (7:automatic, 3:Scotch, 4:PORD, 5:METIS), default: 7
    mmss%mumps_par%icntl(8)  = 7                            ! row and column scaling  7: automatic scaling
    mmss%mumps_par%icntl(18) = 3
    mmss%mumps_par%icntl(14) = 50                           ! MAXS
    
    if (mmss%use_BLR_compression) then
      mmss%mumps_par%icntl(35) = 1                          ! Block-low-rank (BLR) compression. 0: off (default), 1: automatic, 2: factorisation and solution, 3: only factorisation
      mmss%mumps_par%cntl(7)   = mmss%epsilon_BLR                ! Accuracy of BLR approximation
    endif
    
    call DMUMPS(mmss%mumps_par)
      
    mmss%analyzed = .true.
  
  end subroutine mumps_analyze
  
  subroutine mumps_factorize(mmss,a_mat)
  
    use data_structure, only: type_SP_MATRIX
    implicit none
    
    type(type_MUMPS_SOLVER) :: mmss
    type(type_SP_MATRIX)    :: a_mat
    
    mmss%mumps_par%A_loc   => a_mat%val
    mmss%mumps_par%irn_loc => a_mat%irn
    mmss%mumps_par%jcn_loc => a_mat%jcn    
      
    mmss%mumps_par%JOB = 2                                   ! factorisation
    call DMUMPS(mmss%mumps_par)
    
    mmss%analyzed = .true.
    
    return
      
  end subroutine mumps_factorize
  
  subroutine mumps_solve(mmss,rhs_vec)
    use mod_integer_types
    use mpi_mod
    use data_structure, only: type_RHS
    
    implicit none
    
    type(type_MUMPS_SOLVER) :: mmss
    type(type_RHS)          :: rhs_vec
    integer(kind=int_all)   :: k
    integer                 :: ierr
    
    mmss%mumps_par%rhs => rhs_vec%val
    
    mmss%mumps_par%JOB = 3                                   ! Solve

    call DMUMPS(mmss%mumps_par)
    
    call MPI_Bcast(mmss%mumps_par%rhs,mmss%mumps_par%n,MPI_DOUBLE_PRECISION,0,mmss%mumps_par%comm,ierr)    
    
    if (mmss%scaled) then
      do k=1,mmss%mumps_par%n
        rhs_vec%val(k) =  mmss%mumps_par%rhs(k)/mmss%solution_scaling(k)
      enddo
    endif

    return
    
  end subroutine mumps_solve
  
  subroutine mumps_finalize(mmss)
    implicit none
    
    type(type_MUMPS_SOLVER) :: mmss  

    mmss%mumps_par%JOB = -2
    call DMUMPS(mmss%mumps_par)
    
    mmss%comm = 0
    if (associated(mmss%solution_scaling)) deallocate(mmss%solution_scaling)
    mmss%solution_scaling => Null()
    mmss%initialized = .false.
    mmss%analyzed    = .false.
    mmss%equilibrium = .false.    
    mmss%scaled      = .false.
    mmss%refine      = .false.    
    
    return
    
  end subroutine mumps_finalize
  
#endif
end module mod_mumps
