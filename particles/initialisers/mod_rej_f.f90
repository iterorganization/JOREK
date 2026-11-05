!> Module defining abstract interfaces and implemtation for 
!> real space rejection functions used during particle initialisation.
!>
!> A rejection funciton (rej_f) maps field values at a sample point to
!> an acceptance probability in [0,1]. The spatial_pdf derived type bundles 
!> a procedure pointer with the variable list it expects, so both pieces of
!> information always travel together
!>
!> ## Variable index convention (vars(:) entries)
!>
!>  Positive index -> JOREK variable number (model-dependent)
!>  -1             -> R   (major radius)
!>  -2             -> Z   (height)
!>  -3             -> Phi (toroidal angle)
!>
!> The P(:) array passed to a rej_f function is ordered to match vars(:),
!> so a function must be paired with the vars array it was written for. 
!> Use spatial_pdf_from_name() to get a correctly-paired spatial_pdf.

module mod_rej_f
    use equil_info       ! (R_axis,Z_axis), etc..
    use mod_model_settings
    use mpi
    use mod_interp
    use data_structure

    implicit none

    !> Normalisation range for current_pdf, in units of j_tor/R
    !> Set from the equilibrium by calibrate_current_pdf() when the 'current'
    !> pdf is constructed - never use current_pdf without that call.
    real*8 :: jzmin = 0.d0, jzmax = 0.d0

    private

    public :: spatial_pdf
    public :: rej_f
    public :: spatial_pdf_from_name
    public :: check_spatial_pdf
    public :: eval_rej_f
    public :: itpa_tae_pdf, RZ_pdf, analytical_pdf, current_pdf

    ! =============================================================================
    ! Interface
    ! =============================================================================
    
    !> Primary interface for spatial rejection functions used in particle init
    !>
    !> @param n     Number of field values in P and gradP
    !> @param P     Field values at the sample point, ordered as vars(:)
    !> @param gradP (3,n) = (d/dR, d/dZ, d/dphi) of each entry of P.
    !>              Only computed when needs_grad = .true.; zero otherwise.
    !> @return      Acceptance probability in [0,1]
    !>
    !> Must be pure : it is called from inside the OpenMP sampling loops.
    abstract interface
        pure function rej_f(n, P, gradP)
            implicit none
            integer,                intent(in) :: n
            real*8, dimension(n),   intent(in) :: P
            real*8, dimension(3,n), intent(in) :: gradP
            real*4                             :: rej_f
        end function rej_f
    end interface

    ! =============================================================================
    ! spatial_pdf -> Type bundling rejection f and vars it expects
    ! =============================================================================

    !> Bundle type pairing a spatial rejection function with
    !> the JOREK / (R,Z,Phi) variable list it expects
    !>
    !> Always construct via spatial_pdf_from_name() or by setting
    !> both components together - the procedure f is written to 
    !> expect exactly size(vars) values in P(:), ordered to match vars(:).
    !>
    !> Example (custom profile)
    !>
    !>      type(spatial_pdf) :: pdf
    !>      pdf%f          => my_rej_function
    !>      pdf%vars       = [var_psi, -1]   ! psi and R
    !>      pdf%needs_grad = .true.          ! only if my_rej_function requires gradP
    type :: spatial_pdf
        procedure(rej_f), nopass, pointer :: f => null()
        integer,          allocatable     :: vars(:)
        logical                           :: needs_grad = .false. !> Does f require gradients of any vars?
    end type spatial_pdf
    

    contains
    ! =============================================================================
    ! Constructor
    ! =============================================================================

    !> Returns a correctly-paired spatial_pdf type
    !>
    !> Recongnised names:
    !>    'itpa_tae'        - reproduce EP dist in ITPA TAE benchmark
    !>    'RZ'              - weight by 1/(1+r_minor)^2
    !>    'analytical'      - weight by (1-(r/a)^2)^nu
    !>    'current'         - weight by normalised j_tor
    !>    'none'            - no rejection (f => null, vars unallocated)
    !>
    !> node_list/element_list are the equilibrium the pdf is calibrated against.
    !> only 'current' uses them at the moment, but any profile normalised to 
    !> equilibrium fields will need them
    function spatial_pdf_from_name(name, node_list, element_list) result(pdf)
        character(len=*),        intent(in) :: name
        type(type_node_list),    intent(in) :: node_list
        type(type_element_list), intent(in) :: element_list
        type(spatial_pdf)                   :: pdf

        integer :: ierr

        select case(trim(name))

        case('itpa_tae')
            pdf%f          => itpa_tae_pdf
            pdf%vars       =  [var_psi]
            pdf%needs_grad = .false.
        
        case('RZ')
            pdf%f          => RZ_pdf
            pdf%vars       =  [-2, -1]
            pdf%needs_grad = .false.

        case('analytical')
            pdf%f          => analytical_pdf
            pdf%vars       =  [-2, -1]
            pdf%needs_grad = .false.


        case('current')
            !> NOTE : var_zj = 0 in fullMHD (j_tor is not a stored variable)
            !> This rej_f is only valid for models where var_zj > 0 (e.g. model600)
            if (var_zj == 0) then
                write(*,*) "ERROR (mod_rej_f): 'current' spatial PDF requires model"
                write(*,*) "  with var_zj > 0. j_tor is not a stored variable in this model"
                call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
            endif
            call calibrate_current_pdf(node_list, element_list)
            pdf%f          => current_pdf
            pdf%vars       =  [-1, var_zj]
            pdf%needs_grad = .false.


        case('none')
            !> Leave f => null() and vars unallocated
            !> caller should handle this case by skipping use of rej_f
            !> no rej_f should be passed to initialiser, sampling will be uniform in space

        case default
            write(*,*) "ERROR (mod_rej_f): Unknown spatial PDF name '", trim(name), "'"
            write(*,*) "      Valid names: 'itpa_tae', 'RZ', 'analytical', "
            write(*,*) "                   'current', 'none'"
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end select
    end function spatial_pdf_from_name

    !> Set the j_tor/R normalisation range used by current_pdf, once per run.
    !>
    !> find_variable_mimax returns the extrema of one variable over an element's edges
    !> j_tor and R are bounded independently and combined afterwardds, which can only widen the range.
    subroutine calibrate_current_pdf(node_list, element_list)
        use mod_newton_methods, only: find_variable_minmax

        !> i/o vars
        type(type_node_list),    intent(in) :: node_list
        type(type_element_list), intent(in) :: element_list
 
        !> internal vars
        real*8  :: zj_lo, zj_hi, R_lo, R_hi, e_lo, e_hi
        integer :: i_elm, my_id, ierr

        zj_lo = 1.d10; zj_hi = -1.d10
        R_lo  = 1.d10; R_hi  = -1.d10
        
        do i_elm = 1, element_list%n_elements
            call find_variable_minmax(node_list, element_list, i_elm, var_zj, e_lo, e_hi)
            zj_lo = min(zj_lo, e_lo); zj_hi = max(zj_hi, e_hi)
            call find_variable_minmax(node_list, element_list, i_elm, -1, e_lo, e_hi)
            R_lo  = min(R_lo, e_lo);  R_hi  = max(R_hi, e_hi)
        enddo

        !> R>0 on any grid, so j_tor/R is extremal at an R endpoint;
        jzmin = min(zj_lo/R_lo, zj_lo/R_hi)
        jzmax = max(zj_hi/R_lo, zj_hi/R_hi)

        call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

        if (jzmax .le. jzmin) then
            if (my_id == 0) then
                write(*,*) "ERROR (mod_rej_f): 'current' spatial PDF found no spread in j_tor/R"
                write(*,*) " j_tor range ", zj_lo, zj_hi
                write(*,*) " R     range ", R_lo, R_hi
            endif
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif

        if (my_id == 0) then
            write(*,'(A,2ES12.3)') "  j_tor    range : ", zj_lo, zj_hi
            write(*,'(A,2Es12.3)') "  j_tor/R  range : ", jzmin, jzmax
        endif
    end subroutine calibrate_current_pdf

    !> Check that a spatial_pdf satisfies the invariants the samplers rely on.
    !> Does nothing for an inavtive pdf (f => null(), ie init_pdf = 'none')
    subroutine check_spatial_pdf(pdf, caller)
        type(spatial_pdf), intent(in) :: pdf
        character(len=*),  intent(in) :: caller !> name to quote in error message

        integer :: ierr

        if (.not. associated(pdf%f)) return

        if (.not. allocated(pdf%vars)) then
            write(*,*) "ERROR (", caller, "): spatial_pdf has a rejection function but no vars(:)"
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif

        !> Ascendinf order puts the geometric entries (<= 0) first, which is how
        !> the samplers split vars(:) into geometric and JOREK variables
        if (any(pdf%vars(2:) < pdf%vars(:size(pdf%vars)-1))) then
            write(*,*) "ERROR (", caller, "): spatial_pdf vars(:) must be in ascending order, got ", pdf%vars
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif
    end subroutine check_spatial_pdf

    !> Evaluate the variables a spatial rejection function needs at one sample point.
    !>
    !> Fills P(:), and gradP(:,:) when needs_grad, in the order given by variables(:),
    !> following the index convenction documented in mod_rej_f:
    !>     >0: JOREK variable number
    !>      0: constant 1, -1: R, -2: Z, -3: Phi
    !> variables(:) must be sorted ascending, so the n_geom non-positive entries
    !> come first and the n_mhd JOREK variables last.
    !>
    !> gradP is (d/dR, d/dZ, d/dphi) and is set to zero when needs_grad is .false.
    pure function eval_rej_f(node_list, element_list, i_elm, s, t, phi, R, Z, space_pdf) result(f)

    !> i/o vars
    type(type_node_list),    intent(in)  :: node_list
    type(type_element_list), intent(in)  :: element_list
    integer,                 intent(in)  :: i_elm
    real*8,                  intent(in)  :: s, t       !> local element coordinates
    real*8,                  intent(in)  :: R, Z, phi  !> global cylindrical coordinates
    type(spatial_pdf),       intent(in)  :: space_pdf
    real*4                               :: f          !> acceptance probability
    
    !> internal vars
    real*8, parameter :: EPS_JAC = 1.d-12    !> below this element is degenerate
    integer           :: n_geom, n_mhd, k
    real*8            :: P(size(space_pdf%vars)), gradP(3,size(space_pdf%vars))
    real*8            :: P_s(count(space_pdf%vars > 0)), P_t(count(space_pdf%vars > 0))
    real*8            :: P_phi(count(space_pdf%vars > 0))
    real*8            :: R_i, R_S, R_t, Z_i, Z_s, Z_t, xjac, inv_xjac

    n_mhd  = count(space_pdf%vars > 0)
    n_geom = size(space_pdf%vars) - n_mhd
    gradP  = 0.d0

    !> Geomteric vars
    do k = 1, n_geom
        select case (space_pdf%vars(k))
        case(0);  P(k) = 1.d0
        case(-1); P(k) = R   ; if (space_pdf%needs_grad) gradP(1,k) = 1.d0
        case(-2); P(k) = Z   ; if (space_pdf%needs_grad) gradP(2,k) = 1.d0
        case(-3); P(k) = phi ; if (space_pdf%needs_grad) gradP(3,k) = 1.d0
        endselect
    enddo

    !> MHD Variables : values only unless the rej_f asked for gradients
    !> interpolated within the element
    if (n_mhd .ge. 1) then
        if (space_pdf%needs_grad) then
        call interp_PRZ(node_list, element_list, i_elm, space_pdf%vars(n_geom+1:), n_mhd, &
            s, t, phi, P(n_geom+1:), P_s, P_t, P_phi, R_i, R_s, R_t, Z_i, Z_s, Z_t)

        !> interp_PRZ gives s,t derivatives - transform to R,Z
        !> degenerate element (e.g axis point) -> leave gradient zero
        xjac     = R_s*Z_t - R_t*Z_s
        inv_xjac = 0.d0
        if (abs(xjac) > EPS_JAC) inv_xjac = 1.d0/xjac

        do k = 1, n_mhd
            gradP(1,n_geom+k) = ( Z_t*P_s(k) - Z_s*P_t(k)) * inv_xjac
            gradP(2,n_geom+k) = (-R_t*P_s(k) + R_s*P_t(k)) * inv_xjac
            gradP(3,n_geom+k) = P_phi(k)
        enddo
        else
        call interp_PRZ(node_list, element_list, i_elm, space_pdf%vars(n_geom+1:), n_mhd, &
            s, t, phi, P(n_geom+1:), R_i, Z_i)
        endif
    endif

    f = space_pdf%f(size(space_pdf%vars), P, gradP)
    end function eval_rej_f


    ! =============================================================================
    ! Rejection functions
    ! =============================================================================


    !> rejection function to produce EP spatial, distribution from ITPA TAE benchmark
    !> A. Könies et al 2018 Nucl. Fusion 58 126027, https://doi.org/10.1088/1741-4326/aae4e6
    !>
    !> Expected vars : [var_psi]
    !>                 P(1) = psi
    pure function itpa_tae_pdf(n, P, gradP) result(f)
        integer,                intent(in) :: n
        real*8, dimension(n),   intent(in) :: P
        real*8, dimension(3,n), intent(in) :: gradP
        real*4                             :: f
        
        !> internal vars
        real*8              :: s, psi_norm, coeff(0:3)

        coeff(0)=0.49123
        coeff(1)=0.298228
        coeff(2)=0.198739
        coeff(3)=0.521298

        psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)
        s        = 0.957 * psi_norm + 0.043 * psi_norm**2

        f = real(coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2)))), 4)
    end function itpa_tae_pdf

    !> Weight as 1/(1 + r_minor)^2 - broad profile peaked on axis
    !>
    !> Expected vars : [-2, -1]
    !>                 P(1) = Z
    !>                 P(2) = R
    pure function RZ_pdf(n, P, gradP) result(f)
        integer,                intent(in) :: n
        real*8, dimension(n),   intent(in) :: P
        real*8, dimension(3,n), intent(in) :: gradP
        real*4                             :: f

        real*8 :: minor_r

        minor_r = sqrt((P(1) - ES%Z_axis)**2 + (P(2) - ES%R_axis)**2)
        f = real(1.d0 / (1.d0 + minor_r)**2, 4)
    end function RZ_pdf
    
    !> Analytically prescribed profile (1 - (r_minor/a)^2)^nu, nu=2
    !>
    !> Expected vars : [-2, -1]
    !>                 P(1) = Z
    !>                 P(2) = R
    pure function analytical_pdf(n, P, gradP) result(f)
        integer,                intent(in) :: n
        real*8, dimension(n),   intent(in) :: P
        real*8, dimension(3,n), intent(in) :: gradP
        real*4                             :: f

        real*8, parameter :: nu = 2.d0
        real*8            :: minor_r

        minor_r = sqrt((P(1) - ES%Z_axis)**2 + (P(2) - ES%R_axis)**2)
        f = real(max(1.d0 - (minor_r/ES%LCFS_a)**2, 0.d0)**nu, 4)
    end function analytical_pdf

    !> Weight proportional to normalised toroidal current density j_tor
    !>
    !> Expected vars : [-1, var_zj]
    !>                 P(1) = R
    !>                 P(2) = j_tor
    !>
    !> NOTE : Won't work for fMHD models, see note in constructor above
    pure function current_pdf(n, P, gradP) result(f)
        integer,                intent(in) :: n
        real*8, dimension(n),   intent(in) :: P
        real*8, dimension(3,n), intent(in) :: gradP
        real*4                             :: f

        f = real((P(2)/P(1) - jzmin) / (jzmax - jzmin), 4)
        f = min(max(f, 0.0e0), 1.0e0)
    end function current_pdf
end module mod_rej_f
