module mod_sparse_data

#ifdef USE_MUMPS
  use mod_mumps, only:      type_MUMPS_SOLVER
#endif  
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
  use mod_pastix, only:     type_PASTIX_SOLVER
#endif
#ifdef USE_STRUMPACK
  use mod_strumpack, only:  type_STRUMPACK_SOLVER
#endif
  use data_structure, only: type_PRECOND
  
  integer, parameter :: pastix = 1, mumps = 2, strumpack = 3


  type type_SP_SOLVER
#ifdef USE_MUMPS
    type(type_MUMPS_SOLVER)     :: mmss
#endif  
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
    type(type_PASTIX_SOLVER)    :: ptss
#endif
#ifdef USE_STRUMPACK
    type(type_STRUMPACK_SOLVER) :: spss
#endif
    type(type_PRECOND)          :: pc
    
    integer                     :: index_now                           !< current time step index (absolute)
    real(kind=8)                :: tstep                               !< current time step value
    real(kind=8)                :: tstep_prev                          !< previous time step
    integer                     :: istep                               !< curent time step index within jstep group
    
    integer                     :: iter_precon                         !< maximum number of iterations without pc update (input)
    integer                     :: max_steps_noUpdate                  !< maximum number of time steps without pc update (input)
    integer                     :: iter_max                            !< maximum allowed number of iterations (input)
    
    integer                     :: n_since_update = 0                  !< number of time steps since last pc update
    integer                     :: iter_prev = 0                       !< number of iterations in the previous step
    integer                     :: iter_gmres                          !< number of iterations in the current step
    real(kind=8)                :: iter_tol                            !< iterative convergence criteria
    logical                     :: solve_only = .false.                !< flag for updating PC matrix (.true. - no update/factorization needed)
    logical                     :: step_success = .false.              !< flag indicating successfull time step completion
    logical                     :: iterative = .false.                 !< flag indicating use of iterative solver 
    logical                     :: equilibrium = .false.               !< flag indicating equilibrium solver (with duplicate entries in sparse matrix)
    
    integer                     :: library = pastix                    !< solver library (default=pastix)

    logical                     :: verbose = .true.                    !< flag for logfile printout
    
  end type type_SP_SOLVER
   
  private
  public :: type_SP_SOLVER, mumps, pastix, strumpack

end module mod_sparse_data
