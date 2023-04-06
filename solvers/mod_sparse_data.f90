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

  private
  public :: type_SP_SOLVER, mumps, pastix, strumpack

  integer, parameter :: pastix = 1, mumps = 2, strumpack = 3

  type type_newton_solver
    integer                       :: maxNewton    = 20
    real(kind=8)                  :: gamma_Newton = 0.5
    real(kind=8)                  :: alpha_Newton = 2.d0
    real(kind=8), pointer         :: store_value(:,:,:,:) => null()
    real(kind=8), pointer         :: store_delta(:,:,:,:) => null()
    integer                       :: it = 1                   !< current iteration

    contains
    procedure :: newton_setup, newton_finalize
  end type type_newton_solver

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
    type(type_newton_solver)    :: newton

    integer                     :: index_now                           !< current time step index (absolute)
    real(kind=8)                :: tstep                               !< current time step value
    real(kind=8)                :: tstep_prev                          !< previous time step
    integer                     :: istep                               !< curent time step index within jstep group

    integer                     :: iter_precon                         !< maximum number of iterations without pc update (input)
    integer                     :: max_steps_noUpdate                  !< maximum number of time steps without pc update (input)
    integer                     :: iter_max                            !< maximum allowed number of GMRES/BiCGSTAB iterations (input)

    integer                     :: n_since_update = 0                  !< number of time steps since last pc update
    integer                     :: iter_prev = 0                       !< number of iterations in the previous step
    integer                     :: iter_gmres                          !< number of iterations in the current step
    real(kind=8)                :: iter_tol                            !< iterative convergence criteria
    logical                     :: solve_only = .false.                !< flag for updating PC matrix (.true. - no update/factorization needed)
    logical                     :: step_success = .false.              !< flag indicating successfull time step completion
    logical                     :: iterative = .false.                 !< flag indicating use of iterative solver
    logical                     :: equilibrium = .false.               !< flag indicating equilibrium solver (with duplicate entries in sparse matrix)
    logical                     :: use_newton = .false.                !< flag for using iterative Newton method

    integer                     :: library = pastix                    !< solver library (default=pastix)

    logical                     :: verbose = .true.                    !< flag for logfile printout

  contains
    procedure :: setup
    procedure :: finalize

  end type type_SP_SOLVER

  contains

!> Set solver parameters
  subroutine setup(self)
    use phys_module, only: gmres, iter_precon, gmres_max_iter, max_steps_noUpdate, gmres_tol, &
                           use_pastix, use_mumps, use_strumpack, use_newton
    class(type_SP_SOLVER)     :: self

    self%iterative          = gmres
    self%iter_precon        = iter_precon
    self%iter_gmres         = iter_precon
    self%iter_max           = gmres_max_iter
    self%max_steps_noUpdate = max_steps_noUpdate
    self%iter_tol           = gmres_tol
    self%iter_prev          = 0
    self%n_since_update     = 0
    self%use_newton         = use_newton
    if (use_strumpack) then
      self%library = strumpack
    elseif (use_mumps) then
      self%library = mumps
    elseif (use_pastix) then
      self%library = pastix
    endif
    return
  end subroutine setup

!> Finilize solver instance
  subroutine finalize(self)
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
    use mod_pastix, only: pastix_finalize
#endif
#ifdef USE_MUMPS
    use mod_mumps, only: mumps_finalize
#endif
#ifdef USE_STRUMPACK
    use mod_strumpack, only: strumpack_finalize
#endif
    implicit none

    class(type_SP_SOLVER)     :: self

    if (self%verbose) write(*,*) "Finalizing solver"

#if (defined USE_PASTIX) || (defined USE_PASTIX6)
    if (self%ptss%initialized) call pastix_finalize(self%ptss)
#endif
#ifdef USE_MUMPS
    if (self%mmss%initialized) call mumps_finalize(self%mmss)
#endif
#ifdef USE_STRUMPACK
    if (self%spss%initialized) call strumpack_finalize(self%spss)
#endif

    self%solve_only   = .false.
    self%step_success = .false.
    self%iterative    = .false.
    self%equilibrium  = .false.
    self%verbose      = .true.

    call self%newton%newton_finalize()

    return
  end subroutine finalize

!> Set Newton parameters
  subroutine newton_setup(self)
    use phys_module, only: maxNewton, gamma_Newton, alpha_Newton
    class(type_newton_solver)     :: self

    self%maxNewton = maxNewton
    self%gamma_Newton = gamma_Newton
    self%alpha_Newton = alpha_Newton
    return
  end subroutine newton_setup

!> Finalize Newton solver
  subroutine newton_finalize(self)
    class(type_newton_solver)     :: self

    if (associated(self%store_value)) deallocate(self%store_value)
    if (associated(self%store_delta)) deallocate(self%store_delta)

    self%store_value => null()
    self%store_delta => null()
    self%maxNewton = 20
    self%gamma_Newton = 0.5
    self%alpha_Newton = 2.d0
    self%it = 1
    return
  end subroutine newton_finalize

end module mod_sparse_data
